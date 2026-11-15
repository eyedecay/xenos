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
                      half_width_m=300, samples=121):
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
                    "crs": ds.crs.to_string(), "source": ds.tags().get("PRODUCT_ID", "Uploaded lunar DTM")}


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
    return result
