"""Evaluate scan rankings against user-supplied, independent analog-site coordinates."""
from __future__ import annotations
import csv
import io
import math
import statistics


def validate_known_sites(csv_bytes, ranking_pool, match_radius_m=5000.0):
    """Measure nearest sampled-candidate distance and PSI rank for known sites.

    The CSV must have site_name, latitude, longitude columns. A site's optional
    tolerance_m overrides the default match radius. Results are diagnostics, not
    a calibrated accuracy estimate; validation sites must be withheld from tuning.
    """
    if not ranking_pool:
        raise ValueError("Run a scan with candidate-pool diagnostics before validation")
    try:
        text = csv_bytes.decode("utf-8-sig")
        reader = csv.DictReader(io.StringIO(text))
    except Exception as exc:
        raise ValueError(f"Could not read validation CSV: {exc}")
    required = {"site_name", "latitude", "longitude"}
    if not required.issubset(set(reader.fieldnames or [])):
        raise ValueError("Validation CSV needs site_name, latitude, and longitude columns")
    ranked = sorted(ranking_pool,
                    key=lambda item: (item.get("psi") is not None,
                                      item.get("psi") if item.get("psi") is not None
                                      else item.get("terrain_score", -1)), reverse=True)
    rows = []
    for line_number, raw in enumerate(reader, start=2):
        if not any((value or "").strip() for value in raw.values()):
            continue
        try:
            lat, lon = float(raw["latitude"]), float(raw["longitude"])
            tolerance = float(raw.get("tolerance_m") or match_radius_m)
        except (TypeError, ValueError):
            raise ValueError(f"Invalid coordinate or tolerance on CSV row {line_number}")
        if not (-90 <= lat <= 90 and -180 <= lon <= 180 and tolerance > 0):
            raise ValueError(f"Coordinate or tolerance out of range on CSV row {line_number}")
        distances = []
        for item in ranked:
            lat1, lon1 = math.radians(lat), math.radians(lon)
            lat2, lon2 = math.radians(item["lat"]), math.radians(item["lon"])
            dlat, dlon = lat2 - lat1, lon2 - lon1
            a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
            distance = 6_371_000 * 2 * math.asin(min(1.0, math.sqrt(a)))
            distances.append(distance)
        nearest_index = min(range(len(distances)), key=distances.__getitem__)
        nearest = ranked[nearest_index]
        rows.append({
            "site_name": raw["site_name"].strip(), "nearest_candidate_rank": nearest_index + 1,
            "nearest_candidate_region": nearest.get("region"),
            "nearest_candidate_latitude": nearest["lat"],
            "nearest_candidate_longitude": nearest["lon"],
            "distance_m": distances[nearest_index], "tolerance_m": tolerance,
            "within_tolerance": distances[nearest_index] <= tolerance,
            "candidate_psi": nearest.get("psi"),
            "candidate_terrain_score": nearest.get("terrain_score"),
        })
    if not rows:
        raise ValueError("Validation CSV contains no site rows")
    summary = {
        "sites": len(rows),
        "sites_with_candidate_within_tolerance": sum(row["within_tolerance"] for row in rows),
        "fraction_within_tolerance": sum(row["within_tolerance"] for row in rows) / len(rows),
        "median_nearest_distance_m": statistics.median(row["distance_m"] for row in rows),
        "independent_validation_required": "Use sites withheld from target/profile and score tuning",
    }
    return rows, summary
