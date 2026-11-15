"""Regional candidate generation, stratified sampling and ranking."""
from __future__ import annotations
import math
import tempfile
from pathlib import Path
import numpy as np
from ..config import DEM_DIR, GEOLOGY
from .terrain import (circular_zone_masks, riley_terrain_ruggedness, slope_distribution,
                     annular_relief_proxy_score, stratified_candidate_indices)
from .geology import load_geology, geology_composition
from .scoring import (coverage_aware_psi, MIN_GEOLOGY_COVERAGE, MIN_DEM_COVERAGE,
                      MIN_CELLS_ACROSS)

def scan_candidates(diameter_m: int, target_depth_diameter_ratio: float,
                    top_per_tile: int = 72, max_candidates: int = 20,
                    dem_paths=None, geology_path=GEOLOGY,
                    region_name="Columbia River Basalts", relevant_geology_codes=("vb",)):
    """Discover candidates from local DEM sources, stratify terrain sampling, then rank with GLiM."""
    import rasterio
    from pyproj import CRS, Transformer
    from rasterio.warp import transform_bounds
    from scipy.ndimage import convolve

    tile_paths = sorted(Path(path) for path in dem_paths) if dem_paths is not None else sorted(DEM_DIR.glob("*_DEM.tif"))
    if not tile_paths:
        raise ValueError("No DEM GeoTIFFs were provided for this search region")
    bounds_geographic = []
    for path in tile_paths:
        with rasterio.open(path) as ds:
            if ds.crs is None:
                raise ValueError(f"DEM has no CRS: {path.name}")
            bounds_geographic.append(transform_bounds(ds.crs, "EPSG:4326", *ds.bounds))
    center_lon = (min(b[0] for b in bounds_geographic) + max(b[2] for b in bounds_geographic)) / 2
    center_lat = (min(b[1] for b in bounds_geographic) + max(b[3] for b in bounds_geographic)) / 2
    metric_crs = CRS.from_proj4(
        f"+proj=laea +lat_0={center_lat} +lon_0={center_lon} +datum=WGS84 +units=m +no_defs +type=crs"
    )
    geology = load_geology(geology_path, metric_crs)
    to_equal_area = Transformer.from_crs("EPSG:4326", metric_crs, always_xy=True)
    candidates = []
    region_bounds = []
    screened_windows = 0
    for path in tile_paths:
        with rasterio.open(path) as ds:
            region_bounds.append(transform_bounds(ds.crs, "EPSG:4326", *ds.bounds))
            arr = ds.read(1, masked=True).astype("float32").filled(np.nan)
            height, width = arr.shape
            lat = (region_bounds[-1][1] + region_bounds[-1][3]) / 2
            if ds.crs.is_geographic:
                dx = math.hypot(ds.transform.a, ds.transform.d) * 111320 * math.cos(math.radians(lat))
                dy = math.hypot(ds.transform.b, ds.transform.e) * 111320
            else:
                unit_factor = ds.crs.axis_info[0].unit_conversion_factor
                dx = math.hypot(ds.transform.a, ds.transform.d) * unit_factor
                dy = math.hypot(ds.transform.b, ds.transform.e) * unit_factor
            zones = circular_zone_masks(diameter_m, dx, dy)
            masks = {key: value.astype("float32") for key, value in zones.items()}
            valid = np.isfinite(arr)
            surface = np.where(valid, arr, 0).astype("float32")

            def local_mean(values, mask):
                kernel = masks[mask]
                count = convolve(valid.astype("float32"), kernel, mode="constant", cval=0)
                total = convolve(np.where(valid, values, 0).astype("float32"), kernel,
                                 mode="constant", cval=0)
                mean = np.divide(total, count, out=np.full_like(total, np.nan), where=count > 0)
                return mean, count / kernel.sum()

            center_mean, center_coverage = local_mean(surface, "center")
            outer_mean, outer_coverage = local_mean(surface, "outer_annulus")
            _, footprint_coverage = local_mean(surface, "footprint")
            relief = outer_mean - center_mean
            gx, gy = np.gradient(np.where(valid, arr, np.nanmedian(arr)), dy, dx)
            slope = np.degrees(np.arctan(np.hypot(gx, gy))).astype("float32")
            tri = riley_terrain_ruggedness(arr, valid)
            flank_mean_slope, flank_coverage = local_mean(slope, "flank")

            half_rows = zones["footprint"].shape[0] // 2
            half_cols = zones["footprint"].shape[1] // 2
            stride_y = max(2, int(round(half_rows / 2)))
            stride_x = max(2, int(round(half_cols / 2)))
            rows, cols = np.mgrid[half_rows:height - half_rows:stride_y,
                                  half_cols:width - half_cols:stride_x]
            rr, cc = rows.ravel(), cols.ravel()
            ratios = relief[rr, cc] / float(diameter_m)
            keep = (
                np.isfinite(ratios) & (relief[rr, cc] > 0) &
                (footprint_coverage[rr, cc] >= MIN_DEM_COVERAGE) &
                (center_coverage[rr, cc] >= 0.75) &
                (outer_coverage[rr, cc] >= 0.90) &
                (flank_coverage[rr, cc] >= 0.75)
            )
            rr, cc, ratios = rr[keep], cc[keep], ratios[keep]
            screened_windows += len(rr)
            if len(rr) == 0:
                continue
            scores = np.maximum(0.0, 100 * (1 - np.abs(ratios - target_depth_diameter_ratio) /
                                             max(target_depth_diameter_ratio, 0.01)))
            choose = stratified_candidate_indices(scores, top_per_tile)
            for idx in choose:
                row, col = int(rr[idx]), int(cc[idx])
                x, y = rasterio.transform.xy(ds.transform, row, col, offset="center")
                lon, latc = Transformer.from_crs(ds.crs, "EPSG:4326", always_xy=True).transform(x, y)
                east, north = to_equal_area.transform(lon, latc)
                r0, r1 = row - half_rows, row + half_rows + 1
                c0, c1 = col - half_cols, col + half_cols + 1
                plot_window = arr[r0:r1, c0:c1].copy()
                plot_valid = valid[r0:r1, c0:c1]
                plot_window[~(zones["footprint"] & plot_valid)] = np.nan
                local_tri = tri[r0:r1, c0:c1]
                tri_samples = local_tri[zones["footprint"] & plot_valid & np.isfinite(local_tri)]
                local_slope = slope[r0:r1, c0:c1]
                slope_matrix = slope_distribution(local_slope, zones["flank"], plot_valid)
                cells_across = diameter_m / max(dx, dy)
                terrain_score = annular_relief_proxy_score(float(ratios[idx]), target_depth_diameter_ratio)
                candidates.append({
                    "tile": path.stem.replace("Copernicus_DSM_COG_10_", "").replace("_DEM", ""),
                    "region": region_name,
                    "lon": float(lon), "lat": float(latc), "diameter_m": int(diameter_m),
                    "annular_relief_m": float(relief[row, col]),
                    "annular_relief_ratio": float(ratios[idx]),
                    "terrain_score": terrain_score,
                    "flank_mean_slope_deg": float(flank_mean_slope[row, col]),
                    "flank_slope_distribution": slope_matrix,
                    "tri_mean_m": float(tri_samples.mean()) if tri_samples.size else float("nan"),
                    "dem_coverage_fraction": float(footprint_coverage[row, col]),
                    "cells_across": float(cells_across),
                    "pixel_size_x_m": float(dx), "pixel_size_y_m": float(dy),
                    "window": plot_window, "equal_area_x": float(east), "equal_area_y": float(north),
                })

    # Calculate geology across the stratified pool before combined-score ranking.
    for candidate in candidates:
        geo = geology_composition(candidate["lon"], candidate["lat"], diameter_m,
                                  geology, metric_crs=metric_crs,
                                  relevant_codes=relevant_geology_codes)
        candidate["geology"] = geo
        score = coverage_aware_psi(
            candidate["terrain_score"], geo["score"], geo["coverage_fraction"],
            candidate["dem_coverage_fraction"], candidate["cells_across"],
            geo["overlap_fraction"],
        )
        candidate["psi"] = score["psi"]
        candidate["psi_status"] = score["status"]
        candidate["psi_reasons"] = score["reasons"]
    ranking_pool = [{
        "tile": candidate["tile"], "region": candidate["region"],
        "lat": candidate["lat"], "lon": candidate["lon"],
        "terrain_score": candidate["terrain_score"],
        "geology_score": candidate["geology"]["score"],
        "geology_coverage_fraction": candidate["geology"]["coverage_fraction"],
        "dem_coverage_fraction": candidate["dem_coverage_fraction"],
        "psi": candidate["psi"], "psi_eligible": candidate["psi"] is not None,
    } for candidate in candidates]
    candidates.sort(key=lambda c: (c["psi"] is not None,
                                   c["psi"] if c["psi"] is not None else c["terrain_score"]),
                    reverse=True)
    selected = []
    min_separation = diameter_m * 1.5
    for candidate in candidates:
        if all(math.hypot(candidate["equal_area_x"] - prior["equal_area_x"],
                           candidate["equal_area_y"] - prior["equal_area_y"]) >= min_separation
               for prior in selected):
            selected.append(candidate)
            if len(selected) >= max_candidates:
                break
    if region_bounds:
        region_extent = [min(b[0] for b in region_bounds), min(b[1] for b in region_bounds),
                         max(b[2] for b in region_bounds), max(b[3] for b in region_bounds)]
    else:
        region_extent = None
    return selected, {"region_name": region_name,
                      "tiles_read": len(tile_paths),
                      "glim_features": len(geology), "candidate_windows": screened_windows,
                      "candidate_pool_size": len(candidates),
                      "ranking_pool": ranking_pool,
                      "candidate_locations": len(selected),
                      "region_extent_lonlat": region_extent,
                      "min_geology_coverage_for_psi": MIN_GEOLOGY_COVERAGE,
                      "min_dem_coverage_for_psi": MIN_DEM_COVERAGE,
                      "min_cells_across_for_psi": MIN_CELLS_ACROSS}


def scan_uploaded_region(diameter_m, target_depth_diameter_ratio, dem_uploads,
                        geology_name, geology_bytes, region_name="Uploaded Earth region",
                        top_per_tile=72, max_candidates=20,
                        relevant_geology_codes=("vb",)):
    """Scan user-supplied GeoTIFFs and a GeoPackage without modifying source files."""
    import tempfile

    if not dem_uploads:
        raise ValueError("Upload at least one georeferenced DEM GeoTIFF")
    if not geology_bytes or not str(geology_name).lower().endswith((".gpkg", ".geojson")):
        raise ValueError("Upload a geological polygon file (.gpkg or .geojson)")
    with tempfile.TemporaryDirectory(prefix="xenos_region_") as temp_dir:
        temp_root = Path(temp_dir)
        dem_paths = []
        for index, upload in enumerate(dem_uploads):
            name, content = upload
            safe_name = Path(name).name
            if not safe_name.lower().endswith((".tif", ".tiff")):
                raise ValueError(f"DEM must be a GeoTIFF: {safe_name}")
            path = temp_root / f"{index}_{safe_name}"
            path.write_bytes(content)
            dem_paths.append(path)
        geology_path = temp_root / Path(geology_name).name
        geology_path.write_bytes(geology_bytes)
        return scan_candidates(
            diameter_m, target_depth_diameter_ratio, top_per_tile=top_per_tile,
            max_candidates=max_candidates, dem_paths=dem_paths,
            geology_path=geology_path, region_name=region_name,
            relevant_geology_codes=relevant_geology_codes,
        )
