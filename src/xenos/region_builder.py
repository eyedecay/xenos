"""Create additive Earth region packs from public Copernicus tiles and local GLiM."""
from __future__ import annotations

import hashlib
import json
import math
import shutil
import tempfile
from datetime import date
from pathlib import Path

from .config import DATA_DIR

REGION_PACKS_DIR = DATA_DIR / "earth/regions"
GLOBAL_GLIM = DATA_DIR / "earth/geology/glim/LiMW_GIS 2015.gdb"
COPERNICUS_BUCKET = "https://copernicus-dem-30m.s3.amazonaws.com"
GLIM_CITATION = ("Hartmann, J. & Moosdorf, N. (2012). The new global lithological map database GLiM. "
                 "Geochemistry, Geophysics, Geosystems, 13, Q12004. DOI: 10.1029/2012GC004370")


def copernicus_tile_name(latitude_degree: int, longitude_degree: int) -> str:
    """Return the official GLO-30 object prefix for a 1° southwest-corner tile."""
    if not -90 <= latitude_degree <= 89 or not -180 <= longitude_degree <= 179:
        raise ValueError("Copernicus GLO-30 tile origin must be within Earth bounds")
    lat_hemi, lat_value = ("N", latitude_degree) if latitude_degree >= 0 else ("S", -latitude_degree)
    lon_hemi, lon_value = ("E", longitude_degree) if longitude_degree >= 0 else ("W", -longitude_degree)
    return f"Copernicus_DSM_COG_10_{lat_hemi}{lat_value:02d}_00_{lon_hemi}{lon_value:03d}_00_DEM"


def copernicus_tiles_for_bounds(bounds):
    """List all 1° GLO-30 tile IDs intersecting (west, south, east, north)."""
    west, south, east, north = map(float, bounds)
    if not (-180 <= west < east <= 180 and -90 <= south < north <= 90):
        raise ValueError("Use a non-wrapping WGS84 bounding box within lon −180..180, lat −90..90")
    lon_start, lon_stop = math.floor(west), math.ceil(east)
    lat_start, lat_stop = math.floor(south), math.ceil(north)
    return [copernicus_tile_name(lat, lon)
            for lat in range(lat_start, lat_stop)
            for lon in range(lon_start, lon_stop)]


def _clip_global_glim(bounds, output_path, source_path=GLOBAL_GLIM,
                      layer="GLiM_export"):
    """Read only GLiM polygons intersecting the requested area, then clip them."""
    import geopandas as gpd
    import pyogrio
    from pyproj import CRS, Transformer
    from shapely.geometry import box
    from shapely.ops import transform as transform_geometry

    if not Path(source_path).exists():
        raise FileNotFoundError(f"Global GLiM source not found: {source_path}")
    info = pyogrio.read_info(source_path, layer=layer)
    source_crs = CRS.from_user_input(info["crs"])
    west, south, east, north = map(float, bounds)
    wgs_extent = box(west, south, east, north)
    wgs_to_source = Transformer.from_crs("EPSG:4326", source_crs, always_xy=True)
    source_extent = transform_geometry(wgs_to_source.transform, wgs_extent)
    query_bounds = source_extent.bounds
    frame = gpd.read_file(source_path, layer=layer, bbox=query_bounds,
                           columns=["Litho", "xx"], engine="pyogrio")
    if frame.empty:
        raise ValueError("The local GLiM database contains no mapped polygons in this area")
    clip_extent = gpd.GeoSeries([wgs_extent], crs="EPSG:4326").to_crs(frame.crs).iloc[0]
    frame = frame[["Litho", "xx", "geometry"]].copy()
    frame = frame[frame.geometry.notna() & ~frame.geometry.is_empty]
    frame = frame[frame.geometry.intersects(clip_extent)].copy()
    frame.geometry = frame.geometry.intersection(clip_extent)
    frame = frame[frame.geometry.notna() & ~frame.geometry.is_empty]
    if frame.empty:
        raise ValueError("No GLiM polygons remained after clipping to the requested area")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_file(output_path, layer="units", driver="GPKG", engine="pyogrio")
    return len(frame)


def _download_tile(tile_name, destination, timeout=(20, 180)):
    """Download one unsigned public S3 COG. Return False when that tile is absent."""
    import requests

    url = f"{COPERNICUS_BUCKET}/{tile_name}/{tile_name}.tif"
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    try:
        with requests.get(url, stream=True, timeout=timeout) as response:
            if response.status_code == 404:
                return False
            response.raise_for_status()
            with partial.open("wb") as output:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        output.write(chunk)
        if partial.stat().st_size < 1024:
            raise ValueError(f"Downloaded Copernicus tile is unexpectedly small: {tile_name}")
        partial.replace(destination)
    finally:
        partial.unlink(missing_ok=True)
    return True


def _crop_dem_tile(source_path, output_path, bounds):
    """Write an analysis copy cropped to the user's WGS84 search rectangle."""
    import numpy as np
    import rasterio
    from rasterio.mask import mask
    from rasterio.warp import transform_geom
    from shapely.geometry import box, mapping

    with rasterio.open(source_path) as source:
        if source.crs is None:
            raise ValueError(f"Copernicus tile has no CRS: {source_path.name}")
        geometry = transform_geom("EPSG:4326", source.crs,
                                  mapping(box(*map(float, bounds))))
        nodata = -9999.0
        try:
            clipped, transform = mask(source, [geometry], crop=True,
                                      nodata=nodata, filled=True)
        except ValueError as exc:
            raise ValueError(f"The search rectangle did not intersect {source_path.name}") from exc
        data = clipped.astype("float32")
        if not np.isfinite(data).any() or np.all(data == nodata):
            raise ValueError(f"The search rectangle contains no valid elevation in {source_path.name}")
        profile = source.profile.copy()
        profile.update(driver="GTiff", dtype="float32", count=data.shape[0],
                       height=data.shape[1], width=data.shape[2], transform=transform,
                       nodata=nodata, compress="deflate", tiled=False)
        profile.pop("blockxsize", None)
        profile.pop("blockysize", None)
        profile.pop("interleave", None)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(output_path, "w", **profile) as destination:
        destination.write(data)
    return output_path


def build_copernicus_region_pack(name, bounds, *, output_root=REGION_PACKS_DIR,
                                 glim_source=GLOBAL_GLIM, max_tiles=4,
                                 progress_callback=None, place_result=None):
    """Fetch a small Copernicus area and clip local global GLiM into a new pack.

    Source products are read only. The complete pack is staged and renamed into
    place only after all processing succeeds, so failed runs leave no partial pack.
    """
    name = str(name).strip()
    if not name:
        raise ValueError("Enter a name for this Earth region")
    tiles = copernicus_tiles_for_bounds(bounds)
    if len(tiles) > max_tiles:
        raise ValueError(f"This area intersects {len(tiles)} tiles; reduce it to at most {max_tiles} tiles")
    canonical_bounds = [round(float(value), 6) for value in bounds]
    area_hash = hashlib.sha256(json.dumps(canonical_bounds).encode()).hexdigest()[:10]
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    # A changed label must not create another copy of the same geographic data.
    existing_packs = []
    for manifest_path in output_root.glob("*/region.json"):
        try:
            existing = json.loads(manifest_path.read_text())
            old_bounds = [round(float(value), 6) for value in existing.get("bounds_wgs84", [])]
            existing_dems = [path for pattern in existing.get("dem_globs", [])
                             for path in manifest_path.parent.glob(pattern)]
            if (old_bounds == canonical_bounds and existing.get("copernicus_tiles")
                    and existing_dems and all(path.is_file() for path in existing_dems)):
                existing_packs.append((manifest_path, existing))
        except (OSError, ValueError, TypeError):
            continue
    if existing_packs:
        manifest_path, existing = max(
            existing_packs,
            key=lambda item: (bool(item[1].get("place_lookup")),
                              len(str(item[1].get("name", "")))),
        )
        return {"id": existing.get("id", manifest_path.parent.name),
                "name": existing.get("name", name), "path": manifest_path.parent,
                "tiles_downloaded": len(existing.get("copernicus_tiles", [])),
                "tiles_missing": existing.get("unavailable_copernicus_tiles", []),
                "glim_features": int(existing.get("glim_polygon_features_after_clip", 0)),
                "reused": True}
    region_id = f"copernicus-area-{area_hash}"
    final_dir = output_root / region_id
    if final_dir.exists():
        raise FileExistsError(f"Region pack already exists at {final_dir}; select it from Installed region catalog")
    stage = Path(tempfile.mkdtemp(prefix=".xenos-region-", dir=output_root))
    try:
        source_dem_dir = stage / "source_dem"
        analysis_dem_dir = stage / "dem"
        source_dem_dir.mkdir()
        analysis_dem_dir.mkdir()
        available, missing = [], []
        for index, tile in enumerate(tiles, start=1):
            if progress_callback:
                progress_callback(f"Fetching Copernicus GLO-30 tile {index}/{len(tiles)}: {tile}")
            source_path = source_dem_dir / f"{tile}.tif"
            if _download_tile(tile, source_path):
                # Validate that a response is actually a georeferenced raster.
                import rasterio
                with rasterio.open(source_path) as dataset:
                    if dataset.crs is None or dataset.count < 1:
                        raise ValueError(f"Downloaded tile is not a georeferenced DEM: {tile}")
                _crop_dem_tile(source_path, analysis_dem_dir / f"{tile}.tif", bounds)
                available.append(tile)
            else:
                missing.append(tile)
        if not available:
            raise ValueError("Copernicus GLO-30 has no public tile for this area")
        if progress_callback:
            progress_callback("Clipping the local global GLiM polygons to the selected area")
        geology_count = _clip_global_glim(
            bounds, stage / "geology" / "glim_region.gpkg", source_path=glim_source)
        citation = ("Copernicus DEM GLO-30 Public COG, accessed " + date.today().isoformat() +
                    " (https://registry.opendata.aws/copernicus-dem/); " + GLIM_CITATION +
                    "; place lookup via Photon/Komoot using OpenStreetMap data")
        manifest = {
            "id": region_id, "name": name, "dem_globs": ["dem/*.tif"],
            "source_dem_globs": ["source_dem/*.tif"],
            "geology_file": "geology/glim_region.gpkg", "citation": citation,
            "source": "Copernicus GLO-30 public AWS tiles + local global GLiM database",
            "bounds_wgs84": canonical_bounds, "copernicus_tiles": available,
            "unavailable_copernicus_tiles": missing,
            "glim_polygon_features_after_clip": geology_count,
            "place_lookup": place_result,
            "notes": [
                "Analysis DEM files are cropped to the user-selected bounding box; downloaded source tiles are retained separately.",
                "Copernicus GLO-30 is a digital surface model; vegetation and structures may affect elevations.",
                "GLiM is a global-scale lithological map; polygon coverage is not a local field survey.",
                "Only source products intersecting this region were copied or derived; original datasets were not modified.",
            ],
        }
        (stage / "region.json").write_text(json.dumps(manifest, indent=2) + "\n")
        stage.replace(final_dir)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return {"id": region_id, "name": name, "path": final_dir,
            "tiles_downloaded": len(available), "tiles_missing": missing,
            "glim_features": geology_count, "reused": False}
