"""Xenos planetary analog discovery package."""
from .analysis.terrain import (
    riley_terrain_ruggedness, slope_distribution, circular_zone_masks,
    circular_window_metrics, annular_relief_proxy_score, stratified_candidate_indices,
    coarsest_dem_pixel_size_m,
)
from .analysis.geology import (
    GLIM_CLASS_LABELS, load_geology, geology_relevance_score,
    geology_footprint_metrics, geology_composition,
)
from .analysis.scoring import coverage_aware_psi
from .analysis.discovery import scan_candidates, scan_uploaded_region, scan_region_collection
from .references import (
    load_lunar_reference, characterize_lroc_products, lunar_dtm_profile,
    parse_crater_target_profile, planetary_dtm_terrain_metrics,
)
