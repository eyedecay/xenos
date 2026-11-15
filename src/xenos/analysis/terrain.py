"""DEM terrain descriptors and crater-scale screening geometry."""
from __future__ import annotations
import math
import numpy as np

def riley_terrain_ruggedness(elevation, valid=None):
    """Compute Riley TRI: sum of absolute center-to-8-neighbor elevation differences."""
    elevation = np.asarray(elevation, dtype="float32")
    if valid is None:
        valid = np.isfinite(elevation)
    else:
        valid = np.asarray(valid, dtype=bool) & np.isfinite(elevation)
    result = np.zeros_like(elevation, dtype="float32")
    neighbor_count = np.zeros_like(elevation, dtype="uint8")
    padded = np.pad(elevation, 1, mode="edge")
    padded_valid = np.pad(valid, 1, mode="constant", constant_values=False)
    center = elevation
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dx == 0 and dy == 0:
                continue
            neighbor = padded[1 + dy:1 + dy + elevation.shape[0],
                              1 + dx:1 + dx + elevation.shape[1]]
            neighbor_valid = padded_valid[1 + dy:1 + dy + elevation.shape[0],
                                          1 + dx:1 + dx + elevation.shape[1]]
            pair_valid = valid & neighbor_valid
            result += np.where(pair_valid, np.abs(center - neighbor), 0)
            neighbor_count += pair_valid.astype("uint8")
    # Do not score edges or nodata-adjacent cells from a partial neighborhood.
    result[(neighbor_count != 8) | ~valid] = np.nan
    return result


def slope_distribution(slope_degrees, mask, valid=None, bin_width_deg=5):
    """Return a transparent binned slope distribution for a selected terrain zone."""
    values = np.asarray(slope_degrees)
    selected = np.asarray(mask, dtype=bool) & np.isfinite(values)
    if valid is not None:
        selected &= np.asarray(valid, dtype=bool)
    edges = np.arange(0, 90 + bin_width_deg, bin_width_deg, dtype=float)
    counts, edges = np.histogram(values[selected], bins=edges)
    return {"edges_deg": edges.tolist(), "counts": counts.tolist(),
            "sample_count": int(selected.sum()), "bin_width_deg": bin_width_deg}


def circular_zone_masks(diameter_m: float, pixel_x_m: float, pixel_y_m: float):
    """Return metric-circular footprint, center, flank, and outer reference masks."""
    if diameter_m <= 0 or pixel_x_m <= 0 or pixel_y_m <= 0:
        raise ValueError("diameter and pixel dimensions must be positive")
    radius = diameter_m / 2
    half_rows = math.ceil(radius / pixel_y_m)
    half_cols = math.ceil(radius / pixel_x_m)
    yy, xx = np.mgrid[-half_rows:half_rows + 1, -half_cols:half_cols + 1]
    distance = np.hypot(yy * pixel_y_m, xx * pixel_x_m)
    footprint = distance <= radius
    center = distance <= radius * 0.25
    flank = (distance > radius * 0.25) & (distance <= radius * 0.65)
    outer_annulus = (distance >= radius * 0.75) & footprint
    if not center.any() or not flank.any() or not outer_annulus.any():
        raise ValueError("DEM resolution is too coarse to resolve the requested circular zones")
    return {"footprint": footprint, "center": center, "flank": flank,
            "outer_annulus": outer_annulus}


def circular_window_metrics(elevation, valid, pixel_x_m: float, pixel_y_m: float,
                            diameter_m: float, slope_degrees=None):
    """Measure DEM annular relief and optional flank slope inside circular metric zones."""
    zones = circular_zone_masks(diameter_m, pixel_x_m, pixel_y_m)
    if elevation.shape != zones["footprint"].shape or valid.shape != elevation.shape:
        raise ValueError("elevation and validity arrays must match the circular mask dimensions")
    valid = valid.astype(bool) & np.isfinite(elevation)

    def mean_and_coverage(mask, values):
        sample = mask & valid & np.isfinite(values)
        count = int(sample.sum())
        return (float(values[sample].mean()) if count else float("nan"),
                count / int(mask.sum()))

    center_mean, center_coverage = mean_and_coverage(zones["center"], elevation)
    annulus_mean, annulus_coverage = mean_and_coverage(zones["outer_annulus"], elevation)
    _, footprint_coverage = mean_and_coverage(zones["footprint"], elevation)
    flank_slope = None
    flank_coverage = 0.0
    if slope_degrees is not None:
        flank_slope, flank_coverage = mean_and_coverage(zones["flank"], slope_degrees)
    relief = annulus_mean - center_mean if np.isfinite(center_mean) and np.isfinite(annulus_mean) else float("nan")
    return {
        "annular_relief_m": relief,
        "annular_relief_ratio": relief / diameter_m if np.isfinite(relief) else float("nan"),
        "flank_mean_slope_deg": flank_slope,
        "dem_coverage_fraction": footprint_coverage,
        "center_coverage_fraction": center_coverage,
        "outer_annulus_coverage_fraction": annulus_coverage,
        "flank_coverage_fraction": flank_coverage,
        "footprint_mask": zones["footprint"],
    }


def annular_relief_similarity_score(earth_annular_relief_ratio: float, target_annular_relief_ratio: float):
    """Score proportional agreement between like-for-like annular-relief ratios."""
    if not np.isfinite(earth_annular_relief_ratio) or not np.isfinite(target_annular_relief_ratio) or target_annular_relief_ratio <= 0:
        return None
    return float(max(0.0, 100 * (1 - abs(earth_annular_relief_ratio - target_annular_relief_ratio) /
                                      target_annular_relief_ratio)))


def annular_relief_proxy_score(annular_relief_ratio: float, target_depth_diameter_ratio: float):
    """Backward-compatible alias for the ratio similarity heuristic."""
    return annular_relief_similarity_score(annular_relief_ratio, target_depth_diameter_ratio)


def stratified_candidate_indices(scores, limit):
    """Sample the full terrain-score range deterministically instead of keeping only its top."""
    scores = np.asarray(scores, dtype=float)
    if limit <= 0 or not scores.size:
        return np.array([], dtype=int)
    order = np.argsort(scores, kind="stable")
    if len(order) <= limit:
        return order
    selected = np.linspace(0, len(order) - 1, num=limit, dtype=int)
    return order[selected]
