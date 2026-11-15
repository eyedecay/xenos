"""Area-weighted polygon overlays and target-relevant geology scores."""
from __future__ import annotations
from collections import defaultdict
import numpy as np
from .terrain import circular_zone_masks
from ..config import GEOLOGY

GLIM_CLASS_LABELS = {
    "su": "Unconsolidated sediments", "ss": "Siliciclastic sedimentary rocks",
    "py": "Pyroclastics", "sm": "Mixed sedimentary rocks",
    "sc": "Carbonate sedimentary rocks", "ev": "Evaporites",
    "va": "Acid volcanic rocks", "vi": "Intermediate volcanic rocks",
    "vb": "Basic volcanic rocks (usually basalt-type)",
    "pa": "Acid plutonic rocks", "pi": "Intermediate plutonic rocks",
    "pb": "Basic plutonic rocks", "mt": "Metamorphics",
    "wb": "Water bodies", "ig": "Ice and glaciers", "nd": "No data",
}

def load_geology(path=GEOLOGY, target_crs="EPSG:5070"):
    """Read GLiM polygons and project them to an equal-area CRS for overlays."""
    import geopandas as gpd

    try:
        frame = gpd.read_file(path, layer="glim_columbia")
    except Exception:
        frame = gpd.read_file(path)
    missing_columns = {"Litho", "xx", "geometry"} - set(frame.columns)
    if missing_columns:
        raise ValueError("Geology polygons need Litho and xx class-code fields; missing: "
                         + ", ".join(sorted(missing_columns)))
    frame = frame[["Litho", "xx", "geometry"]].dropna(subset=["geometry"])
    return frame.to_crs(target_crs)


def geology_relevance_score(class_shares: dict[str, float], relevant_codes=("vb",)):
    """Score target-relevant mapped classes plus contrast without implying composition identity."""
    codes = tuple(str(code).strip().lower() for code in (relevant_codes or ()) if str(code).strip())
    if not codes:
        return {"score": None, "target_share": None, "contrast_share": None,
                "relevant_codes": [], "highland_anorthosite_supported": False}
    target_share = float(np.clip(sum(class_shares.get(code, 0.0) for code in codes), 0.0, 1.0))
    contrast_share = 1.0 - target_share
    relevance = 100.0 * target_share
    contrast = 200.0 * min(target_share, contrast_share)
    score = 0.5 * relevance + 0.5 * contrast
    return {"score": score, "target_share": target_share,
            "basalt_share": float(np.clip(class_shares.get("vb", 0.0), 0.0, 1.0)),
            "contrast_share": contrast_share, "relevant_codes": list(codes),
            "basalt_relevance_score": relevance, "basalt_contrast_score": contrast,
            "highland_anorthosite_supported": False}


def geology_footprint_metrics(footprint, geology, relevant_codes=("vb",)):
    """Intersect a footprint with GLiM and calculate non-overcounted coverage + class areas."""
    from shapely.ops import unary_union

    try:
        indices = geology.sindex.query(footprint, predicate="intersects")
    except Exception:
        indices = [i for i, geom in enumerate(geology.geometry) if geom.intersects(footprint)]
    detail_areas = defaultdict(float)
    major_areas = defaultdict(float)
    clipped_geometries = []
    for i in indices:
        row = geology.iloc[int(i)]
        intersection = row.geometry.intersection(footprint)
        area = intersection.area
        if area <= 0:
            continue
        clipped_geometries.append(intersection)
        litho = str(row.get("Litho", "nd") or "nd").strip(" _") or "nd"
        major = str(row.get("xx", litho[:2]) or "nd").strip(" _")[:2] or "nd"
        detail_areas[litho] += area
        major_areas[major] += area
    if not clipped_geometries:
        return {"composition": [], "class_shares": {}, "coverage_fraction": 0.0,
                "overlap_fraction": 0.0, "area_m2": 0.0, "score": None,
                "basalt_share": 0.0, "contrast_share": 0.0,
                "highland_anorthosite_supported": False}
    union_area = unary_union(clipped_geometries).area
    summed_area = sum(detail_areas.values())
    footprint_area = footprint.area
    coverage = min(union_area / footprint_area, 1.0) if footprint_area else 0.0
    overlap = max(0.0, summed_area - union_area) / footprint_area if footprint_area else 0.0
    mapped_area = sum(major_areas.values())
    class_shares = {code: area / mapped_area for code, area in major_areas.items()}
    relevance = geology_relevance_score(class_shares, relevant_codes)
    composition = [
        {"litho": code, "class_code": code[:2], "class_label": GLIM_CLASS_LABELS.get(code[:2], "Unclassified"),
         "area_m2": area, "share": area / summed_area}
        for code, area in sorted(detail_areas.items(), key=lambda item: item[1], reverse=True)
    ]
    return {"composition": composition, "class_shares": class_shares,
            "coverage_fraction": coverage, "overlap_fraction": overlap,
            "area_m2": union_area, "score": relevance["score"],
            "basalt_share": relevance.get("basalt_share", 0.0),
            "target_share": relevance.get("target_share"),
            "relevant_codes": relevance.get("relevant_codes", []),
            "contrast_share": relevance.get("contrast_share"),
            "basalt_relevance_score": relevance.get("basalt_relevance_score"),
            "basalt_contrast_score": relevance.get("basalt_contrast_score"),
            "highland_anorthosite_supported": False}


def geology_composition(lon: float, lat: float, diameter_m: float, geology,
                        metric_crs="EPSG:5070", relevant_codes=("vb",)):
    """Calculate area-weighted GLiM composition for a circular metric candidate footprint."""
    from pyproj import Transformer
    from shapely.geometry import Point
    transformer = Transformer.from_crs("EPSG:4326", metric_crs, always_xy=True)
    x, y = transformer.transform(lon, lat)
    footprint = Point(x, y).buffer(diameter_m / 2, quad_segs=32)
    return geology_footprint_metrics(footprint, geology, relevant_codes)
