"""Discovery and validation for locally installed Earth region data packs."""
from __future__ import annotations
import json
from pathlib import Path
from .config import DATA_DIR, DEM_DIR, GEOLOGY

REGION_PACKS_DIR = DATA_DIR / "earth/regions"


def discover_region_profiles():
    """Return only profiles whose real DEM and polygon files are present."""
    columbia_dems = sorted(DEM_DIR.glob("*_DEM.tif")) if DEM_DIR.exists() else []
    profiles = []
    if columbia_dems and GEOLOGY.is_file():
        profiles.append({
            "id": "columbia", "name": "Bundled Columbia River Basalts",
            "dem_paths": columbia_dems, "geology_path": GEOLOGY,
            "relevant_glim_codes": ["vb"], "citation": "Local Copernicus GLO-30 and clipped GLiM inputs",
        })
    if not REGION_PACKS_DIR.exists():
        return profiles
    for manifest_path in sorted(REGION_PACKS_DIR.glob("*/region.json")):
        try:
            raw = json.loads(manifest_path.read_text())
            base = manifest_path.parent.resolve()
            dem_paths = []
            for pattern in raw.get("dem_globs", []):
                dem_paths.extend(sorted(base.glob(pattern)))
            geology_path = (base / raw["geology_file"]).resolve()
            # Do not permit manifests to point outside their own region directory.
            geology_path.relative_to(base)
            if not dem_paths or not geology_path.is_file():
                continue
            if any(not p.is_file() or not p.resolve().is_relative_to(base) for p in dem_paths):
                continue
            codes = raw.get("relevant_glim_codes", [])
            if not isinstance(codes, list):
                continue
            profiles.append({
                "id": raw.get("id", manifest_path.parent.name),
                "name": raw["name"], "dem_paths": dem_paths,
                "geology_path": geology_path, "relevant_glim_codes": codes,
                "citation": raw.get("citation", "Not supplied"),
                "area_bounds": raw.get("bounds_wgs84"),
                "copernicus_tiles": raw.get("copernicus_tiles"),
                "place_lookup": raw.get("place_lookup"),
            })
        except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
            continue
    # Don't let differently named copies of the same Copernicus area over-represent
    # one geographic source in a multi-region scan.
    unique = {}
    for profile in profiles:
        bounds = profile.get("area_bounds")
        tiles = profile.get("copernicus_tiles")
        identity = (("copernicus", tuple(bounds), tuple(tiles))
                    if bounds and tiles else ("name", profile["name"]))
        current = unique.get(identity)
        if current is None:
            profile["duplicate_aliases"] = []
            unique[identity] = profile
            continue
        aliases = current.setdefault("duplicate_aliases", [])
        aliases.extend([current["name"], profile["name"],
                        *profile.get("duplicate_aliases", [])])
        # A geocoder-backed label is stronger provenance; otherwise retain the
        # most informative human label without deleting the source packs.
        current_priority = (bool(current.get("place_lookup")), len(current["name"]))
        profile_priority = (bool(profile.get("place_lookup")), len(profile["name"]))
        if profile_priority > current_priority:
            profile["duplicate_aliases"] = aliases
            unique[identity] = profile
    for profile in unique.values():
        profile["duplicate_aliases"] = sorted(
            set(profile.get("duplicate_aliases", [])) - {profile["name"]}
        )
    return list(unique.values())


def discover_target_profiles():
    """List the bundled McGetchin target and valid cited local Moon/Mars crater profiles."""
    from .references import load_lunar_reference, parse_crater_target_profile

    _baseline, reference = load_lunar_reference()
    morphology = reference["morphology"]
    mcgetchin = {
        "id": "mcgetchin-moon", "name": reference["target"]["name"],
        "profile": {
            "name": reference["target"]["name"], "body": reference["target"]["body"],
            "latitude_deg": reference["target"]["latitude_deg"],
            "longitude_deg": reference["target"]["longitude_deg"],
            "diameter_m": morphology["diameter_m"], "depth_m": morphology["average_depth_m"],
            "depth_diameter_ratio": morphology["average_depth_m"] / morphology["diameter_m"],
            "source_citation": "Bundled McGetchin reference; cross-checked against published study",
            "geological_setting": reference["geological_setting"]["description"],
            "relevant_glim_codes": ["vb"],
            "average_wall_slope_deg": morphology["average_wall_slope_deg"],
            "maximum_wall_slope_deg": morphology["maximum_wall_slope_deg"],
            "average_rim_height_m": morphology["average_rim_height_m"],
        },
    }
    profiles = [mcgetchin]
    profile_dir = DATA_DIR / "targets"
    if not profile_dir.exists():
        return profiles
    for path in sorted(profile_dir.glob("*.json")):
        try:
            profile = parse_crater_target_profile(path.read_bytes())
            profiles.append({"id": path.stem, "name": profile["name"], "profile": profile})
        except (OSError, ValueError):
            continue
    names = set()
    unique = []
    for profile in profiles:
        if profile["name"] not in names:
            names.add(profile["name"])
            unique.append(profile)
    return unique
