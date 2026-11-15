"""Reproducible candidate export serializers."""
from __future__ import annotations
import csv
import io
import json
import math


def candidate_rows(candidates):
    rows = []
    for candidate in candidates:
        geo = candidate.get("geology") or {}
        rows.append({
            "region": candidate.get("region"), "tile": candidate.get("tile"),
            "latitude": candidate.get("lat"), "longitude": candidate.get("lon"),
            "diameter_m": candidate.get("diameter_m"),
            "terrain_similarity_score": candidate.get("terrain_score"),
            "PSI": candidate.get("psi"), "PSI_status": candidate.get("psi_status"),
            "annular_relief_m": candidate.get("annular_relief_m"),
            "annular_relief_ratio": candidate.get("annular_relief_ratio"),
            "DEM_coverage": candidate.get("dem_coverage_fraction"),
            "geology_score_75": geo.get("score"),
            "GLiM_coverage": geo.get("coverage_fraction"),
            "GLiM_overlap": geo.get("overlap_fraction"),
            "GLiM_classes": "; ".join(
                f"{item['litho']}:{item['share']:.4f}" for item in geo.get("composition", [])),
            "cells_across": candidate.get("cells_across"),
            "latitude_longitude_note": "WGS84 point; footprint is circular in regional equal-area analysis",
        })
    return rows


def candidates_csv(candidates):
    rows = candidate_rows(candidates)
    output = io.StringIO()
    fields = list(rows[0]) if rows else ["region", "tile", "latitude", "longitude"]
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def _clean(value):
    if isinstance(value, dict):
        return {key: _clean(item) for key, item in value.items() if key != "window"}
    if isinstance(value, (list, tuple)):
        return [_clean(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def candidates_geojson(candidates):
    features = []
    for row in candidate_rows(candidates):
        feature = {"type": "Feature", "geometry": {
            "type": "Point", "coordinates": [row["longitude"], row["latitude"]]},
            "properties": row}
        features.append(feature)
    return json.dumps(_clean({"type": "FeatureCollection", "features": features}), indent=2,
                      allow_nan=False)


def run_metadata_json(metadata):
    return json.dumps(_clean(metadata), indent=2, sort_keys=True, allow_nan=False, default=str)
