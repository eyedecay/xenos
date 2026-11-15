"""Evaluate scan rankings against user-supplied, independent analog-site coordinates."""
from __future__ import annotations

import csv
import io
import math
import statistics


def _distance_m(lat1, lon1, lat2, lon2):
    """Great-circle distance using a spherical Earth approximation."""
    lat1, lon1, lat2, lon2 = map(math.radians, (lat1, lon1, lat2, lon2))
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return 6_371_000 * 2 * math.asin(min(1.0, math.sqrt(a)))


def validate_known_sites(csv_bytes, ranking_pool, match_radius_m=5000.0, top_k=(1, 5, 10)):
    """Evaluate whether high-ranked samples retrieve held-out analog sites.

    For each known site, this reports the first matching candidate's rank and
    whether any of the top K samples fall within its tolerance. A random-pool
    baseline estimates the expected hit fraction if K samples were selected
    uniformly from the same candidate pool. Because sites/tolerances can
    overlap, these are retrieval diagnostics, not calibrated accuracy.
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
    ranked = sorted((item for item in ranking_pool if item.get("psi") is not None),
                    key=lambda item: item["psi"], reverse=True)
    if not ranked:
        raise ValueError("No candidates passed the PSI coverage gates in this scan")
    rows = []
    for line_number, raw in enumerate(reader, start=2):
        if None in raw:
            raise ValueError(f"Unexpected extra CSV fields on row {line_number}")
        if not any((value or "").strip() for value in raw.values()):
            continue
        try:
            lat, lon = float(raw["latitude"]), float(raw["longitude"])
            tolerance = float(raw.get("tolerance_m") or match_radius_m)
        except (TypeError, ValueError):
            raise ValueError(f"Invalid coordinate or tolerance on CSV row {line_number}")
        if not (-90 <= lat <= 90 and -180 <= lon <= 180 and math.isfinite(tolerance) and tolerance > 0):
            raise ValueError(f"Coordinate or tolerance out of range on CSV row {line_number}")
        distances = [_distance_m(lat, lon, item["lat"], item["lon"]) for item in ranked]
        nearest_index = min(range(len(distances)), key=distances.__getitem__)
        nearest = ranked[nearest_index]
        in_tolerance = [distance <= tolerance for distance in distances]
        first_hit_rank = next((i + 1 for i, hit in enumerate(in_tolerance) if hit), None)
        row = {
            "site_name": raw["site_name"].strip(),
            "nearest_candidate_rank": nearest_index + 1,
            "nearest_candidate_region": nearest.get("region"),
            "nearest_candidate_latitude": nearest["lat"],
            "nearest_candidate_longitude": nearest["lon"],
            "distance_m": distances[nearest_index], "tolerance_m": tolerance,
            "within_tolerance": distances[nearest_index] <= tolerance,
            "first_matching_rank": first_hit_rank,
            "pool_candidate_matches": sum(in_tolerance),
            "candidate_psi": nearest.get("psi"),
            "candidate_terrain_score": nearest.get("terrain_score"),
        }
        for k in top_k:
            if int(k) > 0:
                row[f"hit_at_{int(k)}"] = any(in_tolerance[:int(k)])
        rows.append(row)
    if not rows:
        raise ValueError("Validation CSV contains no site rows")
    summary = {
        "sites": len(rows), "eligible_pool_size": len(ranked),
        "sites_with_candidate_within_tolerance": sum(row["within_tolerance"] for row in rows),
        "fraction_within_tolerance": sum(row["within_tolerance"] for row in rows) / len(rows),
        "median_nearest_distance_m": statistics.median(row["distance_m"] for row in rows),
        "random_pool_hit_fraction": (sum(row["pool_candidate_matches"] / len(ranked) for row in rows)
                                     / len(rows)),
        "random_baseline_note": "Expected hit rate from selecting K candidates uniformly without replacement from this same pool.",
        "independent_validation_required": "Use sites withheld from target/profile and score tuning",
    }
    for k in top_k:
        k = int(k)
        if k <= 0:
            continue
        hits = sum(bool(row.get(f"hit_at_{k}", False)) for row in rows)
        eligible_k = min(k, len(ranked))
        summary[f"hit_at_{k}"] = hits / len(rows)
        n = len(ranked)
        per_site_baselines = []
        for row in rows:
            matches = min(int(row["pool_candidate_matches"]), n)
            misses = n - matches
            no_hit_probability = (math.comb(misses, eligible_k) / math.comb(n, eligible_k)
                                  if misses >= eligible_k and n else 0.0)
            per_site_baselines.append(1 - no_hit_probability)
        summary[f"random_baseline_at_{k}"] = statistics.mean(per_site_baselines)
        baseline = summary[f"random_baseline_at_{k}"]
        summary[f"enrichment_at_{k}"] = (summary[f"hit_at_{k}"] / baseline) if baseline else None
    return rows, summary
