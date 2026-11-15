"""Planetary reference loading and lunar observation processing."""
from __future__ import annotations
import json
import math
import numpy as np
from .config import LUNAR
from .analysis.geology import GLIM_CLASS_LABELS

def load_lunar_reference():
    """Load published target values and provenance supplied with the project."""
    return (json.loads((LUNAR / "crater.json").read_text()),
            json.loads((LUNAR / "lunar_reference.json").read_text()))


def characterize_lroc_products():
    """Read actual PDS EDR image arrays and labels; return stretched thumbnails and metadata."""
    import rasterio

    reference = json.loads((LUNAR / "lunar_reference.json").read_text())
    products = []
    for phase, folder, products_key in [
        ("Pre-impact", "pre-impact", "pre_impact"),
        ("Post-impact", "post-impact", "post_impact"),
    ]:
        for filename in reference["imagery"][products_key]:
            path = LUNAR / "pds" / folder / filename
            with rasterio.open(path) as ds:
                # Render a centered square crop of the narrow swath at preview resolution.
                from rasterio.windows import Window
                crop_rows = min(ds.height, ds.width)
                row0 = max(0, (ds.height - crop_rows) // 2)
                raw = ds.read(1, window=Window(0, row0, ds.width, crop_rows),
                              out_shape=(700, 700), masked=True).astype("float32")
                image = np.asarray(raw.filled(np.nan))
                valid = np.isfinite(image) & (image != 0)
                if valid.any():
                    low, high = np.nanpercentile(image[valid], [2, 98])
                    stretched = np.clip((image - low) / max(high - low, 1), 0, 1)
                    gray = np.nan_to_num(stretched * 255).astype("uint8")
                    thumbnail = np.repeat(gray[:, :, None], 3, axis=2)
                    stats = {
                        "valid_fraction": float(valid.mean()),
                        "p02": float(low),
                        "median": float(np.nanmedian(image[valid])),
                        "p98": float(high),
                        "mean": float(np.nanmean(image[valid])),
                    }
                else:
                    thumbnail = None
                    stats = None
                tags = ds.tags()
                products.append({
                    "phase": phase, "product_id": tags.get("PRODUCT_ID", path.stem),
                    "start_time": tags.get("START_TIME", "Unknown"),
                    "stop_time": tags.get("STOP_TIME", "Unknown"),
                    "instrument": tags.get("INSTRUMENT_NAME", "LROC NAC EDR"),
                    "producer": tags.get("PRODUCER_INSTITUTION_NAME", "Unknown"),
                    "shape": f"{ds.width:,} × {ds.height:,}",
                    "thumbnail": thumbnail, "stats": stats,
                })
    return products


def lunar_dtm_profile(file_bytes, center_lon=67.1765, center_lat=1.3536,
                      half_width_m=300, samples=121, source_name="Uploaded lunar DTM"):
    """Extract a center-crossing profile from an uploaded, mapped lunar DTM GeoTIFF."""
    from pyproj import CRS, Transformer
    from rasterio.io import MemoryFile

    with MemoryFile(file_bytes) as memory_file:
        with memory_file.open() as ds:
            if ds.crs is None:
                raise ValueError("DTM has no map CRS; a georeferenced lunar DTM is required")
            lunar_geographic = CRS.from_proj4(
                "+proj=longlat +a=1737400 +b=1737400 +no_defs +type=crs"
            )
            raster_crs = CRS.from_user_input(ds.crs)
            ellipsoid_radius = raster_crs.ellipsoid.semi_major_metre
            if not math.isclose(ellipsoid_radius, lunar_radius := 1_737_400.0, rel_tol=0.01):
                raise ValueError("DTM CRS is not lunar-centered; use the published lunar NAC DTM")
            if raster_crs.is_geographic:
                to_raster = None
            else:
                to_raster = Transformer.from_crs(lunar_geographic, raster_crs, always_xy=True)
            distances = np.linspace(-half_width_m, half_width_m, samples)
            latitude = center_lat + np.degrees(distances / lunar_radius)
            longitude = np.full_like(latitude, center_lon)
            if to_raster is None:
                xs, ys = longitude, latitude
            else:
                xs, ys = to_raster.transform(longitude.tolist(), latitude.tolist())
            values = []
            for x, y in zip(xs, ys):
                row, col = ds.index(x, y)
                if row < 0 or col < 0 or row >= ds.height or col >= ds.width:
                    values.append(np.nan)
                    continue
                sample = ds.read(1, window=((row, row + 1), (col, col + 1)), masked=True)[0, 0]
                values.append(float(sample) if not np.ma.is_masked(sample) else np.nan)
            profile = np.asarray(values, dtype=float)
            return {"distance_m": distances.tolist(), "elevation_m": profile.tolist(),
                    "valid_fraction": float(np.isfinite(profile).mean()),
                    "pixel_size": [float(abs(ds.transform.a)), float(abs(ds.transform.e))],
                    "crs": ds.crs.to_string(), "source": ds.tags().get("PRODUCT_ID", source_name)}


def parse_crater_target_profile(data):
    """Validate a source-attributed crater profile supplied as JSON."""
    if isinstance(data, (bytes, bytearray)):
        data = json.loads(data.decode("utf-8"))
    elif isinstance(data, str):
        data = json.loads(data)
    if not isinstance(data, dict):
        raise ValueError("Target profile must be a JSON object")
    required = ("name", "body", "diameter_m", "depth_m", "source_citation")
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError("Missing required target fields: " + ", ".join(missing))
    if str(data["body"]).strip().lower() not in {"moon", "mars"}:
        raise ValueError("Target body must be Moon or Mars")
    try:
        diameter = float(data["diameter_m"])
        depth = float(data["depth_m"])
    except (TypeError, ValueError):
        raise ValueError("Target diameter_m and depth_m must be numeric")
    if not np.isfinite(diameter) or not np.isfinite(depth) or diameter <= 0 or depth <= 0 or depth >= diameter:
        raise ValueError("Target measurements require 0 < depth_m < diameter_m")
    if not str(data["source_citation"]).strip():
        raise ValueError("Target profile needs a source citation")
    result = dict(data)
    result["name"] = str(data["name"]).strip()
    result["body"] = str(data["body"]).strip().title()
    result["diameter_m"] = diameter
    result["depth_m"] = depth
    result["depth_diameter_ratio"] = depth / diameter
    codes = data.get("relevant_glim_codes", [])
    if not isinstance(codes, list) or any(not isinstance(code, str) for code in codes):
        raise ValueError("relevant_glim_codes must be a list of GLiM two-letter class codes")
    normalized_codes = [code.strip().lower() for code in codes]
    if any(len(code) != 2 or code not in GLIM_CLASS_LABELS for code in normalized_codes):
        raise ValueError("relevant_glim_codes must use recognized GLiM class codes")
    result["relevant_glim_codes"] = sorted(set(normalized_codes))
    if "latitude_deg" in data or "longitude_deg" in data:
        try:
            latitude, longitude = float(data["latitude_deg"]), float(data["longitude_deg"])
        except (TypeError, ValueError):
            raise ValueError("Target latitude_deg and longitude_deg must be numeric when supplied")
        if not (-90 <= latitude <= 90 and -360 <= longitude <= 360):
            raise ValueError("Target coordinates are outside valid latitude/longitude ranges")
        result["latitude_deg"] = latitude
        result["longitude_deg"] = longitude
    return result


def planetary_dtm_terrain_metrics(file_bytes, diameter_m, center_lon=67.1765,
                                  center_lat=1.3536, body="Moon",
                                  source_name="Uploaded planetary DTM",
                                  comparison_pixel_size_m=None):
    """Measure annular relief from a body-referenced DTM using the Earth footprint zones.

    Requires a projected or geographic raster CRS on the same planetary body. It
    returns None-like coverage metrics rather than inferring values from imagery.
    """
    from pyproj import CRS, Transformer
    from rasterio.io import MemoryFile
    from rasterio.windows import Window, from_bounds
    from rasterio.warp import reproject, Resampling
    from rasterio.transform import from_origin

    radii_m = {"Moon": 1_737_400.0, "Mars": 3_389_500.0}
    body_name = str(body).strip().title()
    if body_name not in radii_m:
        raise ValueError("DTM body must be Moon or Mars")
    if diameter_m <= 0:
        raise ValueError("Footprint diameter must be positive")
    radius_body = radii_m[body_name]
    geographic = CRS.from_proj4(
        f"+proj=longlat +a={radius_body} +b={radius_body} +no_defs +type=crs"
    )
    local = CRS.from_proj4(
        f"+proj=aeqd +lat_0={center_lat} +lon_0={center_lon} "
        f"+a={radius_body} +b={radius_body} +units=m +no_defs +type=crs"
    )
    with MemoryFile(file_bytes) as memfile:
        with memfile.open() as ds:
            if ds.crs is None:
                raise ValueError("DTM has no CRS; use a georeferenced planetary elevation raster")
            src_crs = CRS.from_user_input(ds.crs)
            if not math.isclose(src_crs.ellipsoid.semi_major_metre, radius_body,
                                rel_tol=0.01):
                raise ValueError(f"DTM CRS does not use the {body_name} reference radius")
            source_to_local = Transformer.from_crs(src_crs, local, always_xy=True)
            local_to_source = Transformer.from_crs(local, src_crs, always_xy=True)
            # Transform the square containing the circular footprint to source CRS.
            r = float(diameter_m) / 2
            corner_x = [-r, r, r, -r]
            corner_y = [-r, -r, r, r]
            bx, by = local_to_source.transform(corner_x, corner_y)
            window = from_bounds(min(bx), min(by), max(bx), max(by), ds.transform)
            window = window.round_offsets().round_lengths()
            window = window.intersection(Window(0, 0, ds.width, ds.height))
            if window.width < 1 or window.height < 1:
                raise ValueError("Planetary DTM does not cover the target footprint")
            values = ds.read(1, window=window, masked=True).astype("float64")
            # Estimate native local spacing at the crater center, then compare on
            # no finer than the coarsest selected Earth DEM. Average resampling
            # prevents the 3 m lunar DTM from dominating a ~30 m Earth DEM.
            src_transform = ds.window_transform(window)
            sx, sy = source_to_local.transform(
                [src_transform.c, src_transform.c + src_transform.a,
                 src_transform.c + src_transform.b],
                [src_transform.f, src_transform.f + src_transform.d,
                 src_transform.f + src_transform.e])
            native_spacing = max(math.hypot(sx[1] - sx[0], sy[1] - sy[0]),
                                 math.hypot(sx[2] - sx[0], sy[2] - sy[0]), 1e-6)
            analysis_spacing = max(native_spacing, float(comparison_pixel_size_m or native_spacing))
            cells = max(3, int(math.ceil(float(diameter_m) / analysis_spacing)))
            extent = cells * analysis_spacing / 2
            dst_transform = from_origin(-extent, extent, analysis_spacing, analysis_spacing)
            dst = np.full((cells, cells), np.nan, dtype="float64")
            src_array = np.asarray(values.filled(np.nan), dtype="float64")
            reproject(src_array, dst, src_transform=src_transform, src_crs=ds.crs,
                      src_nodata=np.nan, dst_transform=dst_transform, dst_crs=local,
                      dst_nodata=np.nan, resampling=Resampling.average)
            rows, cols = np.indices(dst.shape, dtype="float64")
            xs = dst_transform.c + (cols + 0.5) * dst_transform.a
            ys = dst_transform.f + (rows + 0.5) * dst_transform.e
            distance = np.hypot(xs, ys)
            valid = np.isfinite(dst)
            footprint = distance <= r
            center = distance <= r * 0.25
            annulus = (distance >= r * 0.75) & footprint
            coverage = lambda zone: float((zone & valid).sum() / zone.sum()) if zone.any() else 0.0
            footprint_coverage = coverage(footprint)
            center_coverage = coverage(center)
            annulus_coverage = coverage(annulus)
            center_samples = dst[center & valid]
            annulus_samples = dst[annulus & valid]
            if not center_samples.size or not annulus_samples.size:
                raise ValueError("DTM has no valid elevations in the center or outer annulus")
            pixel_count_across = float(diameter_m) / analysis_spacing
            relief = float(annulus_samples.mean() - center_samples.mean())
            return {
                "body": body_name, "annular_relief_m": relief,
                "annular_relief_ratio": relief / float(diameter_m),
                "dem_coverage_fraction": footprint_coverage,
                "center_coverage_fraction": center_coverage,
                "outer_annulus_coverage_fraction": annulus_coverage,
                "cells_across": pixel_count_across,
                "native_pixel_size_m": native_spacing,
                "analysis_pixel_size_m": analysis_spacing,
                "source": ds.tags().get("PRODUCT_ID", source_name),
                "crs": ds.crs.to_string(),
                "eligible": (footprint_coverage >= 0.95 and center_coverage >= 0.90
                             and annulus_coverage >= 0.90 and pixel_count_across >= 6
                             and relief > 0),
                "method": "same circular inner disk and outer annulus as Earth candidates",
            }
