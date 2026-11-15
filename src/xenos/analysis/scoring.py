"""Coverage gates and normalized PSI combination."""
from __future__ import annotations
from ..config import (MIN_GEOLOGY_COVERAGE, MIN_DEM_COVERAGE,
                      MIN_CELLS_ACROSS, MAX_GEOLOGY_OVERLAP)

def coverage_aware_psi(terrain_score, geology_score, geology_coverage,
                       dem_coverage, cells_across, overlap_fraction=0.0):
    """Return weighted PSI only when spatial coverage and DEM support pass explicit gates."""
    reasons = []
    if geology_score is None:
        reasons.append("geology score unavailable: target-relevant classes or mapped evidence missing")
    elif geology_coverage < MIN_GEOLOGY_COVERAGE:
        reasons.append(f"GLiM footprint coverage {geology_coverage:.0%} < {MIN_GEOLOGY_COVERAGE:.0%}")
    if dem_coverage < MIN_DEM_COVERAGE:
        reasons.append(f"valid DEM footprint coverage {dem_coverage:.0%} < {MIN_DEM_COVERAGE:.0%}")
    if cells_across < MIN_CELLS_ACROSS:
        reasons.append(f"feature spans {cells_across:.1f} pixels < {MIN_CELLS_ACROSS:.0f}-pixel minimum")
    if overlap_fraction > MAX_GEOLOGY_OVERLAP:
        reasons.append(f"GLiM polygons overlap by {overlap_fraction:.1%}; class areas are ambiguous")
    if reasons:
        return {"psi": None, "status": "withheld", "reasons": reasons}
    flags = []
    if cells_across < 10:
        flags.append(f"coarse DEM support ({cells_across:.1f} cells across footprint)")
    flags.append("GLiM source-map scale is absent from this GeoPackage")
    geology_normalized = geology_score * (100.0 / 75.0)
    psi = 0.65 * terrain_score + 0.35 * geology_normalized
    weights = {"terrain": 0.65, "geology": 0.35}
    return {"psi": float(psi), "weights": weights,
            "status": "screening score; interpret with uncertainty", "reasons": flags}
