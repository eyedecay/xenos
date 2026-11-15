"""Xenos: local planetary analog screening with measured source data."""
from __future__ import annotations

import json
import hashlib
import math
from pathlib import Path

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from src.xenos import (
    characterize_lroc_products, load_lunar_reference,
    lunar_dtm_profile, planetary_dtm_terrain_metrics, parse_crater_target_profile,
    coarsest_dem_pixel_size_m,
    scan_candidates, scan_uploaded_region, scan_region_collection,
)
from src.xenos.terramechanics import bekker_wong_scenario, particle_size_summary

from src.xenos.config import ROOT, LUNAR, MGETCHIN_NAC_DTM
from src.xenos.ui.layout import render_page_chrome, metric, format_coordinates
from src.xenos.ui.provenance import render_provenance
from src.xenos.ui.mobility import render_mobility
from src.xenos.ui.reference_panel import render_reference_panel
from src.xenos.ui.candidate_panel import render_candidate_panel
from src.xenos.ranking import compare_rankings
from src.xenos.ui.validation import render_validation
from src.xenos.ui.exports import render_exports
from src.xenos.ui.briefing import render_candidate_briefing
from src.xenos.regions import discover_region_profiles, discover_target_profiles
from src.xenos.region_builder import (
    build_copernicus_region_pack, copernicus_tiles_for_bounds,
)
from src.xenos.place_search import search_named_places, square_bounds_around_place
BASELINE, REFERENCE = load_lunar_reference()
DEFAULT_MORPH = REFERENCE["morphology"]


@st.cache_data(show_spinner=False)
def cached_lroc_products():
    """Read the large LROC EDR images once, then reuse their summary on reruns."""
    return characterize_lroc_products()


@st.cache_data(ttl=86400, show_spinner=False)
def cached_place_search(query):
    """Cache committed Photon searches so repeated queries do not hit the service."""
    return search_named_places(query)


@st.cache_data(show_spinner=False)
def run_pipeline(diameter_m, target_ratio, region_label, region_name, dem_payloads,
                 geology_payload, relevant_geology_codes, target_metric_name, region_specs):
    if region_label == "Installed region catalog":
        profiles = []
        for profile_id, name, dem_paths, geology_path, citation in region_specs:
            profiles.append({"id": profile_id, "name": name,
                             "dem_paths": [Path(path) for path in dem_paths],
                             "geology_path": Path(geology_path), "citation": citation})
        return scan_region_collection(
            diameter_m, target_ratio, profiles,
            relevant_geology_codes=relevant_geology_codes,
            target_metric_name=target_metric_name,
        )
    geology_name, geology_bytes = geology_payload
    return scan_uploaded_region(diameter_m, target_ratio, dem_payloads,
                                geology_name, geology_bytes, region_name=region_name,
                                relevant_geology_codes=relevant_geology_codes,
                                target_metric_name=target_metric_name)


def main():
    if st.session_state.pop("_xenos_select_installed_once", False):
        st.session_state["region_source_select"] = "Installed region catalog"
    region_to_select = st.session_state.pop("_xenos_region_to_select", None)
    if region_to_select:
        st.session_state["installed_region_selection"] = [region_to_select]
    render_page_chrome(st)
    with st.sidebar:
        st.markdown('<div class="eyebrow">Mission parameters</div>', unsafe_allow_html=True)
        target_profiles = discover_target_profiles()
        profile_ids = [profile["id"] for profile in target_profiles]
        upload_target_id = "__upload_target__"
        target_ids = profile_ids + [upload_target_id]
        target_labels = {profile["id"]: f"{profile['name']} · {profile['profile']['body']}" for profile in target_profiles}
        target_labels[upload_target_id] = "Upload cited crater profile · Moon/Mars"
        target_choice_id = st.selectbox("Planetary reference", target_ids,
                                        format_func=lambda key: target_labels[key])
        target_upload = None
        target_error = None
        is_default_target = target_choice_id == "mcgetchin-moon"
        if target_choice_id == upload_target_id:
            target_upload = st.file_uploader("Cited target profile JSON", type=["json"], key="target_profile")
            try:
                active_target = parse_crater_target_profile(target_upload.getvalue()) if target_upload else None
                if active_target is None:
                    target_error = "Upload a source-cited crater profile to use this target."
            except Exception as exc:
                active_target = None
                target_error = str(exc)
        else:
            active_target = dict(next(profile["profile"] for profile in target_profiles
                                      if profile["id"] == target_choice_id))
        if active_target is None:
            active_target = {
                "name": "No target selected", "body": "Unknown", "diameter_m": DEFAULT_MORPH["diameter_m"],
                "depth_m": DEFAULT_MORPH["average_depth_m"],
                "depth_diameter_ratio": DEFAULT_MORPH["average_depth_m"] / DEFAULT_MORPH["diameter_m"],
                "source_citation": "",
            }
        target_valid = not target_error and (target_choice_id != upload_target_id or target_upload is not None)
        st.markdown(f"### {active_target['name']}")
        st.caption(f"{active_target.get('geological_setting', 'Cited crater profile')} · {active_target['body']} · Source: {active_target['source_citation']}")
        planetary_dtm_upload = st.file_uploader(
            "Optional planetary DTM GeoTIFF for like-for-like terrain comparison",
            type=["tif", "tiff"], key="lunar_dtm",
            help="Requires a body-centered CRS and known crater-center coordinates. Xenos measures the same circular annular relief on this DTM and Earth DEM.",
        )
        if planetary_dtm_upload is not None:
            planetary_dtm_bytes = planetary_dtm_upload.getvalue()
            planetary_dtm_source = planetary_dtm_upload.name
        elif is_default_target and MGETCHIN_NAC_DTM.is_file():
            planetary_dtm_bytes = MGETCHIN_NAC_DTM.read_bytes()
            planetary_dtm_source = "NASA LROC NAC 3 m DTM, McGetchin (Zenodo 21481275)"
            st.caption("NASA LROC NAC McGetchin DTM found locally; it will be measured as the target terrain reference.")
        else:
            planetary_dtm_bytes = None
            planetary_dtm_source = None
        st.markdown("---")
        st.markdown("**Earth search regions**")
        available_regions = discover_region_profiles()
        region_source = st.selectbox(
            "Region data",
            ["Installed region catalog", "Explore an area with Copernicus AWS", "Upload one region for this run"],
            key="region_source_select",
        )
        region_profiles = []
        dem_uploads = geology_upload = None
        custom_region_name = "Uploaded Earth region"
        custom_region_citation = ""
        prepare_copernicus_region = False
        if region_source == "Installed region catalog":
            region_names = [profile["name"] for profile in available_regions]
            selected_region_names = st.multiselect(
                "Search one or more installed regions", region_names,
                default=region_names[:1],
                help="Each region is analyzed against its own DEM and geology layer; results are merged with per-region representation.",
                key="installed_region_selection",
            )
            region_profiles = [profile for profile in available_regions if profile["name"] in selected_region_names]
            duplicate_names = sorted({alias for profile in region_profiles
                                       for alias in profile.get("duplicate_aliases", [])})
            if duplicate_names:
                st.info("Some saved entries share the same Copernicus area and are counted once: "
                        + ", ".join(duplicate_names))
            region_choice = "Installed region catalog"
            if available_regions:
                installed = ", ".join(profile["name"] for profile in available_regions)
                st.caption(f"Installed region packs: {installed}. Scores use measured terrain and mapped geology.")
            else:
                st.info("No complete local region packs found. Add a region pack or upload one for this run.")
        elif region_source == "Upload one region for this run":
            custom_region_name = st.text_input("Region name", value="Uploaded Earth region")
            custom_region_citation = st.text_input("DEM and geology source citation", placeholder="Dataset title, version, URL, license")
            dem_uploads = st.file_uploader("Upload georeferenced DEM GeoTIFF tiles", type=["tif", "tiff"], accept_multiple_files=True, key="region_dems")
            geology_upload = st.file_uploader("Upload GLiM-compatible geology polygons", type=["gpkg", "geojson"], key="region_geology", help="Vector layer must carry Litho and xx class-code fields.")
            region_choice = "Uploaded region"
        else:
            st.caption("Search place names, select a match, and Xenos will center a bounded search area there. Place names come from Photon using OpenStreetMap data; Copernicus itself provides tile IDs, not named places.")
            with st.form("place_lookup_form", clear_on_submit=False):
                place_query = st.text_input(
                    "Search for a place", placeholder="Atacama Desert, Haughton Crater, Iceland…",
                    help="Submit a place name to see matching locations. Your query is sent to Photon/Komoot and cached for one day.",
                )
                search_places = st.form_submit_button("Find matching places", width="stretch")
            if search_places:
                try:
                    found_places = cached_place_search(place_query)
                    st.session_state["xenos_place_search_results"] = found_places
                    st.session_state["xenos_place_search_query"] = place_query
                    st.session_state.pop("xenos_place_search_error", None)
                except Exception as exc:
                    st.session_state["xenos_place_search_results"] = []
                    st.session_state["xenos_place_search_error"] = str(exc)
            if st.session_state.get("xenos_place_search_error"):
                st.warning(f"Place search failed: {st.session_state['xenos_place_search_error']}")
            place_results = st.session_state.get("xenos_place_search_results", [])
            if place_results:
                place_options = {item["id"]: item for item in place_results}
                result_set_key = hashlib.sha1("|".join(place_options).encode()).hexdigest()[:10]
                selected_place_id = st.selectbox(
                    "Matching places", list(place_options),
                    format_func=lambda key: place_options[key]["label"],
                    key=f"xenos_place_match_{result_set_key}",
                )
                selected_place = place_options[selected_place_id]
                area_width_km = st.selectbox(
                    "Search area width around this place", [25, 50, 100], index=1,
                    format_func=lambda value: f"{value} km × {value} km",
                    key=f"xenos_area_width_{selected_place_id}",
                    help="This is a square centered on the selected place, not the full boundary of the named feature.",
                )
                copernicus_bounds = square_bounds_around_place(
                    selected_place["longitude"], selected_place["latitude"], area_width_km)
                copernicus_region_name = st.text_input(
                    "Region display name", value=selected_place["name"],
                    key=f"xenos_region_label_{selected_place_id}",
                    help="This label does not change the searched data.",
                )
                selected_tiles = copernicus_tiles_for_bounds(copernicus_bounds)
                st.caption(
                    f"Center: {selected_place['latitude']:.4f}, {selected_place['longitude']:.4f} · "
                    f"area bounds: W {copernicus_bounds[0]:.4f}, S {copernicus_bounds[1]:.4f}, "
                    f"E {copernicus_bounds[2]:.4f}, N {copernicus_bounds[3]:.4f} · "
                    f"Copernicus tiles: {len(selected_tiles)}"
                )
                st.caption("Place search: Photon by Komoot · © OpenStreetMap contributors. Place queries are sent to the geocoding service when you press Find matching places.")
                if len(selected_tiles) > 4:
                    st.warning(f"This area intersects {len(selected_tiles)} Copernicus tiles. Choose a smaller search width (maximum four tiles).")
                with st.form("copernicus_prepare_form", clear_on_submit=False):
                    prepare_copernicus_region = st.form_submit_button(
                        "Download DEM + prepare geology", width="stretch",
                        disabled=len(selected_tiles) > 4,
                    )
                place_result_metadata = {
                    "provider": "Photon/Komoot with OpenStreetMap data",
                    "display_name": selected_place["label"],
                    "latitude": selected_place["latitude"],
                    "longitude": selected_place["longitude"],
                    "osm_type": selected_place["osm_type"],
                    "osm_id": selected_place["osm_id"],
                }
            else:
                copernicus_region_name = ""
                copernicus_bounds = None
                place_result_metadata = None
                st.info("Search for a place to choose an area. You no longer need to enter longitude and latitude manually.")
            region_choice = "Copernicus AWS region"
        st.markdown("**Candidate footprint**")
        if target_error:
            st.error(target_error)
        with st.form("candidate_search_form", clear_on_submit=False):
            diameter = st.slider("Diameter (m)", 100, 5000, int(min(max(active_target["diameter_m"], 100), 5000)), 10)
            st.caption("Searches crater morphology and mapped geology for field-exercise candidates; not an operational site assessment.")
            st.caption("Uses matching annular-relief metrics when a suitable planetary DTM is available; otherwise uses a clearly labeled unlike-measure proxy.")
            run_scan = st.form_submit_button("◉  Run regional analog scan", use_container_width=True)
        st.markdown('<p class="smallnote">Candidate windows are Earth analog candidates, not confirmed craters or landing-site safety assessments.</p>', unsafe_allow_html=True)

    if prepare_copernicus_region:
        try:
            with st.status("Preparing the new Earth region…", expanded=True) as prep_status:
                result = build_copernicus_region_pack(
                    copernicus_region_name, copernicus_bounds,
                    progress_callback=lambda message: prep_status.update(label=message, state="running", expanded=True),
                    place_result=place_result_metadata,
                )
                prep_status.update(
                    label=f"Region ready · {result['tiles_downloaded']} DEM tiles · {result['glim_features']:,} GLiM polygons",
                    state="complete", expanded=False,
                )
            st.session_state["_xenos_select_installed_once"] = True
            st.session_state["_xenos_region_to_select"] = result["name"]
            st.session_state["_xenos_region_created"] = result
            st.rerun()
        except Exception as exc:
            st.error(f"Could not prepare the Copernicus region: {exc}")

    created_region = st.session_state.pop("_xenos_region_created", None)
    if created_region:
        missing_note = (f" · unavailable DEM tiles: {', '.join(created_region['tiles_missing'])}"
                        if created_region["tiles_missing"] else "")
        action = "Reusing the existing" if created_region.get("reused") else "Added"
        st.success(f"{action} {created_region['name']} region: {created_region['tiles_downloaded']} Copernicus DEM tile(s), {created_region['glim_features']:,} clipped GLiM polygons{missing_note}. Review the source limitations, then run the candidate scan.")

    dem_payloads = [(file.name, file.getvalue()) for file in (dem_uploads or [])]
    geology_payload = ((geology_upload.name, geology_upload.getvalue())
                       if geology_upload is not None else ("", b""))
    selected_dem_paths = [path for profile in region_profiles for path in profile["dem_paths"]]
    try:
        earth_comparison_resolution_m = coarsest_dem_pixel_size_m(
            selected_dem_paths if region_source == "Installed region catalog" else (),
            dem_payloads if region_source != "Installed region catalog" else (),
        )
    except Exception as exc:
        earth_comparison_resolution_m = None
        st.sidebar.warning(f"Could not determine Earth DEM resolution for target comparison: {exc}")

    # The planetary DTM only replaces the unlike-measure proxy when the same annular
    # geometry has adequate coverage and resolution on the planetary raster.
    target_ratio = active_target["depth_diameter_ratio"]
    target_metric_name = "published depth/diameter ratio (unlike-measure proxy)"
    target_dtm_metrics = None
    target_dtm_profile = None
    target_dtm_error = None
    if planetary_dtm_bytes is not None:
        try:
            dtm_bytes = planetary_dtm_bytes
            target_lon = active_target.get("longitude_deg", REFERENCE["target"]["longitude_deg"] if is_default_target else None)
            target_lat = active_target.get("latitude_deg", REFERENCE["target"]["latitude_deg"] if is_default_target else None)
            if target_lon is None or target_lat is None:
                raise ValueError("Add latitude_deg and longitude_deg to the cited target profile before uploading its DTM")
            target_dtm_metrics = planetary_dtm_terrain_metrics(
                dtm_bytes, diameter, target_lon, target_lat, active_target["body"],
                source_name=planetary_dtm_source,
                comparison_pixel_size_m=earth_comparison_resolution_m)
            if target_dtm_metrics["eligible"]:
                target_ratio = target_dtm_metrics["annular_relief_ratio"]
                target_metric_name = "resolution-matched DTM annular-relief ratio"
            else:
                target_dtm_error = ("DTM coverage, relief sign, or resolution does not pass the reference gates; "
                                    "the scan uses the published depth/diameter proxy.")
            if active_target["body"] == "Moon":
                target_dtm_profile = lunar_dtm_profile(
                    dtm_bytes, target_lon, target_lat, source_name=planetary_dtm_source)
        except Exception as exc:
            target_dtm_error = str(exc)
            target_dtm_metrics = None
            target_dtm_profile = None
    if target_dtm_error:
        st.sidebar.warning(target_dtm_error)

    MORPH = {
        "diameter_m": active_target["diameter_m"], "average_depth_m": active_target["depth_m"],
        "average_wall_slope_deg": active_target.get("average_wall_slope_deg"),
        "maximum_wall_slope_deg": active_target.get("maximum_wall_slope_deg"),
        "average_rim_height_m": active_target.get("average_rim_height_m"),
    }

    a, b, c, d = st.columns(4)
    with a: metric(st, "Reference diameter", f"{active_target['diameter_m']:g}" if target_valid else "—", "m")
    with b: metric(st, "Target comparison ratio", f"{target_ratio:.3f}" if target_valid else "—", "DTM" if target_metric_name.startswith("DTM") else "published proxy")
    with c: metric(st, "Reference wall slope median", str(active_target.get('average_wall_slope_deg') if active_target.get('average_wall_slope_deg') is not None else 'n/a'), "°")
    with d: metric(st, "Reference imagery", "LROC NAC" if is_default_target else "Not bundled", "")

    region_specs = tuple(
        (p["id"], p["name"], tuple(str(path) for path in p["dem_paths"]),
         str(p["geology_path"]), p.get("citation", "Not supplied"))
        for p in region_profiles
    )
    region_ready = (bool(region_profiles) if region_source == "Installed region catalog"
                    else bool(dem_payloads and geology_payload[1])
                    if region_source == "Upload one region for this run" else False)
    requested_scan_args = (diameter, target_ratio, region_choice,
                           tuple((profile_id, name, dem_paths, geology_path,
                                  tuple((Path(path).stat().st_size, Path(path).stat().st_mtime_ns) for path in dem_paths),
                                  Path(geology_path).stat().st_size, Path(geology_path).stat().st_mtime_ns)
                                 for profile_id, name, dem_paths, geology_path, _citation in region_specs),
                           custom_region_name, custom_region_citation,
                           tuple((name, hashlib.sha256(content).hexdigest()) for name, content in dem_payloads),
                           (geology_payload[0], hashlib.sha256(geology_payload[1]).hexdigest()),
                           target_choice_id, active_target["name"], active_target["body"], active_target["source_citation"],
                           tuple(active_target.get("relevant_glim_codes", [])), target_metric_name,
                           hashlib.sha256(planetary_dtm_bytes).hexdigest() if planetary_dtm_bytes else None)
    scan_ready = region_ready and target_valid
    st.session_state.setdefault("candidates", [])
    st.session_state.setdefault("scan_summary", None)
    st.session_state.setdefault("scan_error", None)
    st.session_state.setdefault("scan_args", None)
    if run_scan and scan_ready:
        scan_status = None
        try:
            with st.status("Reading DEM windows and intersecting GLiM geology…", expanded=True) as scan_status:
                found, scan_summary = run_pipeline(diameter, target_ratio, region_choice, custom_region_name,
                                                   dem_payloads, geology_payload,
                                                   tuple(active_target.get("relevant_glim_codes", [])),
                                                   target_metric_name, region_specs)
                scan_status.update(label=f"Scan complete · retained {len(found)} candidate locations",
                                   state="complete", expanded=False)
            st.session_state.candidates = found
            st.session_state.scan_summary = scan_summary
            st.session_state.scan_error = None
            st.session_state.scan_args = requested_scan_args
        except Exception as exc:
            if scan_status is not None:
                scan_status.update(label="Regional analog scan failed", state="error", expanded=True)
            st.session_state.candidates = []
            st.session_state.scan_summary = None
            st.session_state.scan_error = str(exc)
            st.session_state.scan_args = requested_scan_args
    elif not target_valid:
        st.info("Upload a valid, source-cited crater target profile before scanning.")
    elif not region_ready:
        st.info("Select at least one complete installed region, or upload a georeferenced DEM and geology polygon file.")
    elif st.session_state.get("scan_args") is None:
        st.info("Choose the target, region, and footprint size, then run the scan to discover candidates.")
    candidates_default_state = st.session_state.get("candidates", [])
    candidates = candidates_default_state if scan_ready else []
    for source_index, item in enumerate(candidates):
        item["_source_index"] = source_index
    filtered_candidates = candidates
    sort_key = lambda item: item["psi"] if item["psi"] is not None else -1
    filtered_candidates.sort(key=sort_key, reverse=True)
    saved_scan_args = st.session_state.get("scan_args")
    results_are_stale = saved_scan_args is not None and saved_scan_args != requested_scan_args
    if results_are_stale:
        st.warning("The displayed candidates use the previous scan settings. Run the scan again to apply the current sidebar values.")

    left, right = st.columns([1.02, 1.3], gap="medium")
    render_reference_panel(st, left, is_default_target, REFERENCE, MORPH, LUNAR, go,
                           cached_lroc_products, active_target, target_dtm_profile,
                           target_dtm_metrics, target_metric_name)

    active_candidate = render_candidate_panel(st, right, go, np, math, filtered_candidates, format_coordinates)
    render_candidate_briefing(st, active_candidate)
    if active_candidate:
        st.caption("Potential use: comparative field exercises in crater-like landform mapping and mapped geological context. Xenos does not assess field access, safety, lighting, weather, communications, resources, or rover mobility.")
    st.markdown('<div class="panel"><div class="panel-title">Explainable similarity and mapped composition</div><div class="panel-sub">Terrain and geology components remain visible as independent scores</div>', unsafe_allow_html=True)
    if filtered_candidates:
        candidate = active_candidate
        tcol, gcol = st.columns([.8, 1.2], gap="large")
        with tcol:
            st.metric("Terrain similarity score", f"{candidate['terrain_score']:.0f}/100")
            st.caption(f"Earth annular-relief ratio: {candidate['annular_relief_ratio']:.3f} · target comparison ratio: {target_ratio:.3f} ({target_metric_name})")
            with st.expander("Additional terrain measurements", expanded=False):
                st.metric("Riley terrain ruggedness", f"{candidate['tri_mean_m']:.1f} m")
                st.caption(f"Mean 8-neighbor TRI inside the circular footprint · DEM cells {candidate['cells_across']:.1f} across · valid coverage {candidate['dem_coverage_fraction']:.1%}.")
            geo_score = candidate["geology"]["score"]
            geo_coverage = candidate["geology"]["coverage_fraction"]
            relevance_label = "GLiM basalt relevance + contrast" if is_default_target else "Target-class relevance + contrast"
            geo_label = relevance_label if geo_coverage >= 0.90 else "Provisional geology score · partial coverage"
            geo_value = "—" if geo_score is None else f"{geo_score:.0f}/75"
            st.metric(geo_label, geo_value)
            normalized_geo = None if geo_score is None else geo_score * (100 / 75)
            norm_text = 'unavailable' if normalized_geo is None else f'{normalized_geo:.0f}/100 for PSI'
            if is_default_target:
                st.caption(f"`vb` basic-volcanic share: {candidate['geology']['basalt_share']:.0%} · mapped contrast share: {candidate['geology']['contrast_share']:.0%} · normalized geology component {norm_text}. No GLiM class is treated as anorthosite.")
            else:
                codes = ", ".join(active_target.get("relevant_glim_codes", [])) or "none configured"
                target_share = candidate["geology"].get("target_share")
                share_text = "unavailable" if target_share is None else f"{target_share:.0%}"
                st.caption(f"Target-relevant GLiM codes: {codes} · mapped share {share_text} · normalized geology component {norm_text}. Supply evidence-based codes in the cited target profile.")
            with st.expander("Slope details", expanded=False):
                st.metric("Earth flank-band mean slope", f"{candidate['flank_mean_slope_deg']:.1f}°")
                st.caption("Descriptive only; not scored against the target's median interior wall slope, which uses a different statistic and spatial baseline.")
            st.metric("Planetary similarity index", "Unavailable" if candidate["psi"] is None else f"{candidate['psi']:.0f}/100")
            weights_note = " Formula: 65% terrain similarity + 35% normalized geology." if candidate["psi"] is not None else "."
            st.caption(candidate["psi_status"] + weights_note)
            if candidate["psi_reasons"]:
                st.caption(" · ".join(candidate["psi_reasons"]))
        with gcol:
            composition = candidate["geology"]["composition"]
            if composition:
                frame = pd.DataFrame(composition)
                frame["share_pct"] = frame["share"] * 100
                chart = px.bar(frame.head(8), x="share_pct", y="litho", orientation="h", text=frame.head(8)["share_pct"].map(lambda v: f"{v:.1f}%"), color="share_pct", color_continuous_scale=["#315b50", "#c5f36b"])
                chart.update_layout(height=270, margin=dict(l=8,r=8,t=6,b=6), paper_bgcolor="#11191d", plot_bgcolor="#11191d", font_color="#9eaaa9", xaxis_title="Area share within candidate footprint", yaxis_title="GLiM Litho code", coloraxis_showscale=False)
                chart.update_yaxes(autorange="reversed")
                st.plotly_chart(chart, use_container_width=True, config={"displayModeBar":False})
                st.caption(f"Area-weighted polygon intersections · {candidate['geology']['area_m2']:,.0f} m² mapped · GLiM areal coverage {candidate['geology']['coverage_fraction']*100:.1f}% · polygon overlap {candidate['geology']['overlap_fraction']*100:.2f}%.")
                st.dataframe(pd.DataFrame([{"GLiM Litho code": r["litho"], "Major class": f"{r['class_code']} · {r['class_label']}", "Intersected area (m²)": round(r["area_m2"], 1), "Area share (%)": round(r["share"]*100, 2)} for r in composition]), hide_index=True, use_container_width=True)
                st.caption("GLiM source-map scale is not included in this GeoPackage. Full polygon coverage is not a statement of map accuracy at crater scale.")
            else:
                st.info("No GLiM polygon intersected this candidate footprint; no geology score was assigned.")
                st.caption(f"GLiM areal coverage: {candidate['geology']['coverage_fraction']:.1%}.")
            histogram = candidate["flank_slope_distribution"]
            edges = histogram["edges_deg"]
            centers = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
            with st.expander("Slope distribution and comparison limits", expanded=False):
                slope_fig = go.Figure(go.Bar(x=centers, y=histogram["counts"], marker_color="#6e9b83"))
                if is_default_target:
                    slope_fig.add_vrect(x0=24, x1=41, fillcolor="#c5f36b", opacity=0.12, line_width=0,
                                        annotation_text="McGetchin median-to-max context", annotation_position="top left")
                slope_fig.update_layout(height=210, margin=dict(l=8,r=8,t=20,b=8), paper_bgcolor="#11191d", plot_bgcolor="#11191d", font_color="#9eaaa9", xaxis_title="Earth flank-band slope (°)", yaxis_title="Pixel count")
                st.plotly_chart(slope_fig, use_container_width=True, config={"displayModeBar":False})
                slope_context = "McGetchin 24–41° range is context only" if is_default_target else "No slope range supplied for this custom target"
                st.caption(f"Slope distribution at {histogram['bin_width_deg']}° bins; {histogram['sample_count']} pixels. {slope_context}; not an equivalent filter or score.")
        if is_default_target:
            geology_recipe = "Geology score = 0.5 × (100 × `vb` share) + 0.5 × (200 × min(`vb` share, non-`vb` share)); maximum 75. `vb` is GLiM basic volcanic, usually basalt-type. Other classes are generic mapped contrast; none is assumed to represent anorthosite."
        elif active_target.get("relevant_glim_codes"):
            codes_text = ", ".join(active_target["relevant_glim_codes"])
            geology_recipe = f"Geology score uses target-relevant GLiM codes ({codes_text}) as the relevant share and all other mapped classes as contrast: 0.5 × (100 × relevant share) + 0.5 × (200 × min(relevant share, contrast share)); maximum 75."
        else:
            geology_recipe = "Geology composition is reported, but no target-relevant GLiM classes were supplied; the geology score and PSI are withheld."
        with st.expander("Score formulas and scientific limits", expanded=False):
            st.markdown(f"""**Scoring recipe**
    Terrain similarity score = max(0, 100 × (1 − |Earth annular-relief ratio − target comparison ratio| / target comparison ratio)). The target comparison ratio comes from the same circular annular-relief calculation on a local or uploaded planetary DTM when available; otherwise it uses published depth/diameter as a unlike-measure screening proxy. Annular relief is the mean elevation in the outer 0.75–1.00 radius ring minus the mean elevation in the inner 0–0.25 radius disk. It is **not** measured crater depth. The Earth flank-band slope is descriptive and excluded because it is not methodologically comparable to the target median wall slope.
    {geology_recipe}
    PSI first normalizes the 0–75 GLiM geology score to 0–100. PSI = 0.65 × terrain proxy + 0.35 × normalized geology. PSI is withheld if GLiM coverage is below 90%, valid DEM coverage below 95%, the footprint spans fewer than 6 DEM pixels, or overlapping GLiM polygons exceed 1%. These are screening gates, not a calibrated uncertainty model.""")
    else:
        st.info("Run the regional candidate scan to inspect terrain and GLiM scores.")
    st.markdown('</div>', unsafe_allow_html=True)

    with st.expander("Advanced: candidate diagnostics and validation", expanded=False):
        if filtered_candidates:
            compare_options = list(range(len(filtered_candidates)))
            compare_selected = st.multiselect(
                "Compare candidates (up to 5)", compare_options,
                default=compare_options[:min(3, len(compare_options))], max_selections=5,
                format_func=lambda i: f"{filtered_candidates[i]['tile']} · {format_coordinates(filtered_candidates[i]['lat'], filtered_candidates[i]['lon'])}",
            )
            if compare_selected:
                compare_rows = []
                for i in compare_selected:
                    item = filtered_candidates[i]
                    compare_rows.append({
                        "Candidate": f"{item['tile']} · {format_coordinates(item['lat'], item['lon'], 4)}",
                        "PSI": item["psi"], "Terrain /100": item["terrain_score"],
                        "Geology /75": item["geology"]["score"],
                        "DEM coverage": item["dem_coverage_fraction"],
                        "GLiM coverage": item["geology"]["coverage_fraction"],
                        "Cells across": item["cells_across"],
                    })
                st.dataframe(pd.DataFrame(compare_rows), hide_index=True, use_container_width=True)
            terrain_weight = st.slider("Terrain weight for sensitivity check", 0.0, 1.0, 0.65, 0.05,
                                       help="Exploratory weight variation; does not overwrite the baseline PSI.")
            ranking_pool = (st.session_state.scan_summary or {}).get("ranking_pool", [])
            stability, rankings = compare_rankings(ranking_pool, terrain_weight)
            if rankings:
                st.markdown("**How many of the top ten stay on the shortlist?**")
                st.dataframe(pd.DataFrame(stability), hide_index=True, use_container_width=True)
                st.caption("All scenarios use the same PSI-eligible pool before map thinning. A tie at the cutoff means pool order can change shortlist membership.")
                with st.expander("Inspect scenario ranks"):
                    rows = [{"Scenario": label, "Rank": rank, "Region": item["region"],
                             "Candidate": format_coordinates(item["lat"], item["lon"], 4),
                             "Score": score} for label, rank, item, score in rankings]
                    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
            else:
                st.info("No samples passed the PSI evidence gates for a weight comparison.")
        else:
            st.info("Ranking diagnostics become available after candidates pass the current coverage filter.")
        validation_pool = (st.session_state.scan_summary or {}).get("ranking_pool", [])
        scanned = (st.session_state.scan_summary or {}).get("region_summaries", [])
        odessa_scan = (len(scanned) == 1 and scanned[0].get("region_name") == "Odessa crater study area"
                       and bool(scanned[0].get("region_extent_lonlat")))
        if odessa_scan:
            west, south, east, north = scanned[0]["region_extent_lonlat"]
            odessa_scan = west <= -102.4791667 <= east and south <= 31.7561111 <= north
        if not results_are_stale:
            render_validation(st, validation_pool, odessa_scan)

    render_provenance(st)

    region_citations = ([{"region": p["name"], "citation": p.get("citation", "Not supplied")} for p in region_profiles]
                        if region_source == "Installed region catalog" else
                        [{"region": custom_region_name, "citation": custom_region_citation or "Not supplied by uploader"}])
    run_manifest = {
        "application": "Xenos planetary analog discovery",
        "target": {"name": active_target["name"], "body": active_target["body"],
                   "source_citation": active_target.get("source_citation"),
                   "diameter_m": active_target.get("diameter_m"),
                   "depth_m": active_target.get("depth_m"),
                   "latitude_deg": active_target.get("latitude_deg", REFERENCE["target"]["latitude_deg"] if is_default_target else None),
                   "longitude_deg": active_target.get("longitude_deg", REFERENCE["target"]["longitude_deg"] if is_default_target else None),
                   "relevant_glim_codes": active_target.get("relevant_glim_codes", [])},
        "terrain_comparison": {"ratio": target_ratio, "method": target_metric_name,
                                "dtm_metrics": target_dtm_metrics},
        "search": {"diameter_m": diameter, "regions": region_citations,
                   "sorted_by": "PSI when eligible; terrain score when PSI is withheld",
                   "formula": "65% terrain similarity + 35% normalized geology",
                   "psi_gates": {"minimum_glim_coverage": 0.90, "minimum_dem_coverage": 0.95,
                                 "minimum_cells_across": 6, "maximum_polygon_overlap": 0.01}},
        "input_signatures_and_scan_parameters": requested_scan_args,
        "scan_summary": st.session_state.get("scan_summary"),
        "limitations": ["PSI is a screening index, not a probability or validated analog determination.",
                        "Source rasters and per-pixel arrays are not embedded in the export."],
    }
    render_exports(st, filtered_candidates, run_manifest)

    with st.expander("Experimental: soil and rover scenario calculator", expanded=False):
        render_mobility(st, pd, np, go, particle_size_summary, bekker_wong_scenario)
