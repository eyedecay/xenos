"""Sentinel-2 L2A chip alignment and exploratory spectral-contact metrics."""
from __future__ import annotations

import math
import re

import numpy as np


REQUIRED_BANDS = ("B02", "B03", "B04", "B11", "B12")


def spectral_indices(bands):
    """Return visible iron-oxide proxy and SWIR contrast from aligned reflectance bands."""
    blue, red, swir1, swir2 = (np.asarray(bands[key], dtype="float32")
                               for key in ("B02", "B04", "B11", "B12"))
    iron_proxy = (red - blue) / np.maximum(np.abs(red + blue), 1e-6)
    swir_contrast = (swir1 - swir2) / np.maximum(np.abs(swir1 + swir2), 1e-6)
    return iron_proxy, swir_contrast


def contact_gradient_score(iron_proxy, swir_contrast, valid, threshold_per_pixel=0.01):
    """Score the share of valid pixels with a strong joint spectral-index gradient."""
    iron_proxy = np.asarray(iron_proxy, dtype="float32")
    swir_contrast = np.asarray(swir_contrast, dtype="float32")
    valid = np.asarray(valid, dtype=bool) & np.isfinite(iron_proxy) & np.isfinite(swir_contrast)
    if threshold_per_pixel <= 0:
        raise ValueError("spectral contact threshold must be positive")
    if min(iron_proxy.shape) < 2:
        return {"score": None, "gradient": np.zeros_like(iron_proxy),
                "contact_mask": np.zeros_like(valid), "sample_count": int(valid.sum())}
    fe_y, fe_x = np.gradient(iron_proxy)
    swir_y, swir_x = np.gradient(swir_contrast)
    gradient = np.sqrt(fe_y**2 + fe_x**2 + swir_y**2 + swir_x**2)
    contact = valid & (gradient >= threshold_per_pixel)
    count = int(valid.sum())
    score = 100.0 * float(contact.sum()) / count if count else None
    return {"score": score, "gradient": gradient, "contact_mask": contact,
            "sample_count": count, "threshold_per_pixel": threshold_per_pixel}


def process_sentinel2_uploads(uploaded_files, reflectance_scale=0.0001):
    """Align uploaded L2A B02/B03/B04/B11/B12 GeoTIFF chips onto the B02 grid."""
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.io import MemoryFile
    from rasterio.warp import reproject

    if reflectance_scale <= 0:
        raise ValueError("reflectance scale must be positive")
    file_by_band = {}
    for item in uploaded_files:
        name = getattr(item, "name", "")
        match = re.search(r"(?:^|[_-])(B\d{2}|SCL)(?:[_\.-]|$)", name.upper())
        if match:
            file_by_band[match.group(1)] = item.getvalue()
    missing = sorted(set(REQUIRED_BANDS) - set(file_by_band))
    if missing:
        raise ValueError("Missing Sentinel-2 bands: " + ", ".join(missing))

    reference_memory = MemoryFile(file_by_band["B02"])
    reference = reference_memory.open()
    if reference.crs is None:
        reference.close()
        reference_memory.close()
        raise ValueError("B02 GeoTIFF needs a projected CRS and geotransform")
    if max(reference.width, reference.height) > 2048:
        reference.close()
        reference_memory.close()
        raise ValueError("Upload a cropped Sentinel-2 chip no larger than 2048 pixels per side")
    if not reference.crs.is_projected:
        reference.close()
        reference_memory.close()
        raise ValueError("Upload Sentinel-2 bands in a projected CRS, such as the tile UTM zone")
    shape = (reference.height, reference.width)
    aligned = {}
    valid_masks = []
    file_mems = [reference_memory]
    try:
        for band in REQUIRED_BANDS:
            memory = reference_memory if band == "B02" else MemoryFile(file_by_band[band])
            if band != "B02":
                file_mems.append(memory)
            with memory.open() as source:
                if source.crs is None:
                    raise ValueError(f"{band} has no CRS")
                raw = source.read(1, masked=True).astype("float32").filled(np.nan)
                nodata = source.nodata
                raw_valid = np.isfinite(raw)
                if nodata is not None:
                    raw_valid &= raw != nodata
                source_values = np.where(raw_valid, raw * reflectance_scale, np.nan).astype("float32")
                if band == "B02":
                    values = source_values
                    valid = raw_valid
                else:
                    values = np.full(shape, np.nan, dtype="float32")
                    reproject(source_values, values,
                              src_transform=source.transform, src_crs=source.crs,
                              src_nodata=np.nan, dst_transform=reference.transform,
                              dst_crs=reference.crs, dst_nodata=np.nan,
                              resampling=Resampling.bilinear)
                    valid = np.isfinite(values)
                aligned[band] = values
                valid_masks.append(valid)
        scl_metadata = "not supplied"
        if "SCL" in file_by_band:
            memory = MemoryFile(file_by_band["SCL"])
            file_mems.append(memory)
            with memory.open() as source:
                if source.crs is None:
                    raise ValueError("SCL band has no CRS")
                scl = np.zeros(shape, dtype="uint8")
                reproject(source.read(1), scl,
                          src_transform=source.transform, src_crs=source.crs,
                          dst_transform=reference.transform, dst_crs=reference.crs,
                          resampling=Resampling.nearest)
                # Exclude no-data, saturated pixels, cloud shadow, clouds, cirrus and snow.
                valid_masks.append(~np.isin(scl, [0, 1, 3, 7, 8, 9, 10, 11]))
                scl_metadata = "applied (cloud/shadow/snow classes masked)"
        valid = np.logical_and.reduce(valid_masks)
        indices = spectral_indices(aligned)
        true_color_reflectance = np.stack([aligned["B04"], aligned["B03"], aligned["B02"]], axis=-1)
        rgb = np.zeros(true_color_reflectance.shape, dtype="uint8")
        for channel in range(3):
            channel_values = true_color_reflectance[:, :, channel]
            samples = channel_values[valid & np.isfinite(channel_values)]
            if samples.size:
                low, high = np.percentile(samples, [2, 98])
                stretched = np.clip((channel_values - low) / max(high - low, 1e-6), 0, 1)
                rgb[:, :, channel] = np.nan_to_num(stretched * 255).astype("uint8")
        return {
            "bands": aligned, "valid": valid, "iron_proxy": indices[0],
            "swir_contrast": indices[1], "true_color": rgb,
            "transform": reference.transform, "crs": reference.crs.to_string(),
            "shape": shape, "scl": scl_metadata,
            "scale": reflectance_scale,
        }
    finally:
        reference.close()
        for memory in set(file_mems):
            memory.close()


def sentinel_candidate_metrics(processed, lon, lat, diameter_m, threshold_per_pixel=0.01):
    """Measure Sentinel coverage/contact heterogeneity over a DEM candidate's circle."""
    import rasterio
    from pyproj import CRS, Transformer
    from rasterio.transform import xy

    crs = CRS.from_user_input(processed["crs"])
    transformer = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    x, y = transformer.transform(lon, lat)
    transform = processed["transform"]
    center_row, center_col = rasterio.transform.rowcol(transform, x, y)
    unit_factor = crs.axis_info[0].unit_conversion_factor
    pixel_x = math.hypot(transform.a, transform.d) * unit_factor
    pixel_y = math.hypot(transform.b, transform.e) * unit_factor
    radius = diameter_m / 2
    half_rows = math.ceil(radius / pixel_y)
    half_cols = math.ceil(radius / pixel_x)
    r0, r1 = center_row - half_rows, center_row + half_rows + 1
    c0, c1 = center_col - half_cols, center_col + half_cols + 1
    h, w = processed["shape"]
    clipped = (max(r0, 0), min(r1, h), max(c0, 0), min(c1, w))
    cr0, cr1, cc0, cc1 = clipped
    if cr0 >= cr1 or cc0 >= cc1:
        return {"score": None, "coverage_fraction": 0.0, "contact_mask": None,
                "iron_proxy": None, "swir_contrast": None}
    rows, cols = np.mgrid[cr0:cr1, cc0:cc1]
    xs, ys = xy(transform, rows, cols, offset="center")
    xs = np.asarray(xs).reshape(rows.shape)
    ys = np.asarray(ys).reshape(rows.shape)
    distance = np.hypot((xs - x) * unit_factor, (ys - y) * unit_factor)
    circle = distance <= radius
    local_valid = processed["valid"][cr0:cr1, cc0:cc1] & circle
    coverage = float(local_valid.sum() / circle.sum()) if circle.any() else 0.0
    iron = processed["iron_proxy"][cr0:cr1, cc0:cc1]
    swir = processed["swir_contrast"][cr0:cr1, cc0:cc1]
    signal = contact_gradient_score(iron, swir, local_valid, threshold_per_pixel)
    return {"score": signal["score"] if coverage >= 0.90 else None,
            "coverage_fraction": coverage, "contact_mask": signal["contact_mask"],
            "iron_proxy": iron, "swir_contrast": swir,
            "threshold_per_pixel": threshold_per_pixel,
            "sample_count": signal["sample_count"],
            "window_bounds": (cr0, cr1, cc0, cc1)}
