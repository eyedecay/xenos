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
    characterize_lroc_products, coverage_aware_psi, load_lunar_reference,
    lunar_dtm_profile, parse_crater_target_profile, scan_candidates,
    scan_uploaded_region,
)
from src.xenos.sentinel import process_sentinel2_uploads, sentinel_candidate_metrics
from src.xenos.terramechanics import bekker_wong_scenario, particle_size_summary

from src.xenos.config import ROOT, LUNAR
from src.xenos.ui.layout import render_page_chrome, metric, format_coordinates
from src.xenos.ui.provenance import render_provenance
from src.xenos.ui.mobility import render_mobility
from src.xenos.ui.reference_panel import render_reference_panel
from src.xenos.ui.candidate_panel import render_candidate_panel
BASELINE, REFERENCE = load_lunar_reference()
DEFAULT_MORPH = REFERENCE["morphology"]

render_page_chrome(st)
with st.sidebar:
    st.markdown('<div class="eyebrow">Mission parameters</div>', unsafe_allow_html=True)
    target_choice = st.selectbox("Planetary reference", ["McGetchin Crater · Moon", "Upload cited crater profile · Moon/Mars"])
    target_upload = None
    target_error = None
    if target_choice.startswith("Upload"):
        target_upload = st.file_uploader("Cited target profile JSON", type=["json"], key="target_profile")
        try:
            active_target = parse_crater_target_profile(target_upload.getvalue()) if target_upload else None
            if active_target is None:
                target_error = "Upload a source-cited crater profile to use this target."
        except Exception as exc:
            active_target = None
            target_error = str(exc)
    else:
        active_target = {
            "name": REFERENCE["target"]["name"], "body": REFERENCE["target"]["body"],
            "diameter_m": DEFAULT_MORPH["diameter_m"], "depth_m": DEFAULT_MORPH["average_depth_m"],
            "depth_diameter_ratio": DEFAULT_MORPH["average_depth_m"] / DEFAULT_MORPH["diameter_m"],
            "source_citation": "Bundled McGetchin reference; cross-checked with published study",
            "geological_setting": REFERENCE["geological_setting"]["description"],
            "relevant_glim_codes": ["vb"],
            "average_wall_slope_deg": DEFAULT_MORPH["average_wall_slope_deg"],
            "maximum_wall_slope_deg": DEFAULT_MORPH["maximum_wall_slope_deg"],
            "average_rim_height_m": DEFAULT_MORPH["average_rim_height_m"],
        }
    if active_target is None:
        active_target = {
            "name": "No target selected", "body": "Unknown", "diameter_m": DEFAULT_MORPH["diameter_m"],
            "depth_m": DEFAULT_MORPH["average_depth_m"],
            "depth_diameter_ratio": DEFAULT_MORPH["average_depth_m"] / DEFAULT_MORPH["diameter_m"],
            "source_citation": "",
        }
    is_default_target = target_choice.startswith("McGetchin")
    target_valid = not target_error and (is_default_target or target_upload is not None)
    target_ratio = active_target["depth_diameter_ratio"]
    st.markdown(f"### {active_target['name']}")
    st.caption(f"{active_target.get('geological_setting', 'Cited crater profile')} · {active_target['body']} · Source: {active_target['source_citation']}")
    st.markdown("---")
    st.markdown("**Earth search region**")
    region_choice = st.selectbox("Region", ["Bundled Columbia River Basalts", "Upload another Earth region"])
    dem_uploads = geology_upload = None
    custom_region_name = "Uploaded Earth region"
    if region_choice.startswith("Upload"):
        custom_region_name = st.text_input("Region name", value="Uploaded Earth region")
        dem_uploads = st.file_uploader("Upload georeferenced DEM GeoTIFF tiles", type=["tif", "tiff"], accept_multiple_files=True, key="region_dems")
        geology_upload = st.file_uploader("Upload GLiM-compatible geology polygons", type=["gpkg", "geojson"], key="region_geology", help="Vector layer must carry Litho and xx class-code fields.")
    st.markdown("**Candidate footprint**")
    diameter = st.slider("Diameter (m)", 100, 5000, int(min(max(active_target["diameter_m"], 100), 5000)), 10)
    st.caption("For crater targets, the Earth annular-relief ratio is compared with the target's cited depth/diameter ratio; they remain different physical measurements.")
    st.markdown("---")
    if target_error:
        st.error(target_error)
    run_scan = st.button("◉  Run regional analog scan", use_container_width=True)
    minimum_coverage = st.slider("Minimum DEM and GLiM coverage", 0.0, 1.0, 0.0, 0.05, help="Filters the currently discovered candidate list.")
    sort_choice = st.selectbox("Rank candidates by", ["PSI", "Terrain score", "Geology score"])
    st.markdown('<p class="smallnote">Crater-profile mode screens positive annular-relief patterns. User-supplied DEM and geology layers enable additional Earth regions; outputs remain candidate analogs, not confirmed equivalents.</p>', unsafe_allow_html=True)

MORPH = {
    "diameter_m": active_target["diameter_m"], "average_depth_m": active_target["depth_m"],
    "average_wall_slope_deg": active_target.get("average_wall_slope_deg"),
    "maximum_wall_slope_deg": active_target.get("maximum_wall_slope_deg"),
    "average_rim_height_m": active_target.get("average_rim_height_m"),
}

a, b, c, d = st.columns(4)
with a: metric(st, "Reference diameter", f"{active_target['diameter_m']:g}" if target_valid else "—", "m")
with b: metric(st, "Reference depth / diameter", f"{target_ratio:.2f}" if target_valid else "—", "published")
with c: metric(st, "Reference wall slope median", str(active_target.get('average_wall_slope_deg') if active_target.get('average_wall_slope_deg') is not None else 'n/a'), "°")
with d: metric(st, "Reference imagery", "LROC NAC" if is_default_target else "Not bundled", "")

@st.cache_data(show_spinner=False)
def run_pipeline(diameter_m, depth_diameter_ratio, region_label, region_name, dem_payloads, geology_payload, relevant_geology_codes):
    if region_label == "Bundled Columbia River Basalts":
        return scan_candidates(diameter_m, depth_diameter_ratio, region_name="Columbia River Basalts",
                               relevant_geology_codes=relevant_geology_codes)
    geology_name, geology_bytes = geology_payload
    return scan_uploaded_region(diameter_m, depth_diameter_ratio, dem_payloads,
                                geology_name, geology_bytes, region_name=region_name,
                                relevant_geology_codes=relevant_geology_codes)

dem_payloads = [(file.name, file.getvalue()) for file in (dem_uploads or [])]
geology_payload = ((geology_upload.name, geology_upload.getvalue())
                   if geology_upload is not None else ("", b""))
region_ready = region_choice.startswith("Bundled") or bool(dem_payloads and geology_payload[1])
requested_scan_args = (diameter, target_ratio, region_choice,
                       custom_region_name,
                       tuple((name, hashlib.sha256(content).hexdigest()) for name, content in dem_payloads),
                       (geology_payload[0], hashlib.sha256(geology_payload[1]).hexdigest()),
                       active_target["name"], active_target["body"], active_target["source_citation"],
                       tuple(active_target.get("relevant_glim_codes", [])))
scan_ready = region_ready and target_valid
st.session_state.setdefault("candidates", [])
st.session_state.setdefault("scan_summary", None)
st.session_state.setdefault("scan_error", None)
st.session_state.setdefault("scan_args", None)
if (run_scan or st.session_state.get("scan_args") != requested_scan_args) and scan_ready:
    try:
        with st.spinner("Reading DEM windows and intersecting GLiM geology…"):
            found, scan_summary = run_pipeline(diameter, target_ratio, region_choice, custom_region_name,
                                               dem_payloads, geology_payload,
                                               tuple(active_target.get("relevant_glim_codes", [])))
        st.session_state.candidates = found
        st.session_state.scan_summary = scan_summary
        st.session_state.scan_error = None
        st.session_state.scan_args = requested_scan_args
    except Exception as exc:
        st.session_state.candidates = []
        st.session_state.scan_summary = None
        st.session_state.scan_error = str(exc)
        st.session_state.scan_args = requested_scan_args
elif not target_valid:
    st.info("Upload a valid, source-cited crater target profile before scanning.")
elif not region_ready:
    st.info("Upload at least one georeferenced DEM and one geology polygon file to scan this region.")
candidates_default_state = st.session_state.get("candidates", [])
candidates = candidates_default_state if scan_ready else []
for source_index, item in enumerate(candidates):
    item["_source_index"] = source_index
filtered_candidates = [item for item in candidates
                       if item["dem_coverage_fraction"] >= minimum_coverage
                       and item["geology"]["coverage_fraction"] >= minimum_coverage]
sort_key = {"PSI": lambda item: item["psi"] if item["psi"] is not None else -1,
            "Terrain score": lambda item: item["terrain_score"],
            "Geology score": lambda item: item["geology"]["score"] if item["geology"]["score"] is not None else -1}[sort_choice]
filtered_candidates.sort(key=sort_key, reverse=True)
results_are_stale = st.session_state.get("scan_args") != requested_scan_args
if results_are_stale:
    st.warning("The displayed candidates use the previous scan settings. Run the scan again to apply the current sidebar values.")

left, right = st.columns([1.02, 1.3], gap="medium")
render_reference_panel(st, left, is_default_target, REFERENCE, MORPH, LUNAR, go, lunar_dtm_profile, characterize_lroc_products, active_target)

sentinel_metric = None
candidate, sentinel_metric = render_candidate_panel(st, right, go, np, math, filtered_candidates, candidates, format_coordinates, process_sentinel2_uploads, sentinel_candidate_metrics)


active_candidate = dict(candidate) if filtered_candidates else None
if filtered_candidates and sentinel_metric is not None:
    psi_result = coverage_aware_psi(
        candidate["terrain_score"], candidate["geology"]["score"],
        candidate["geology"]["coverage_fraction"], candidate["dem_coverage_fraction"],
        candidate["cells_across"], candidate["geology"]["overlap_fraction"],
        sentinel_score=sentinel_metric["score"],
        sentinel_coverage=sentinel_metric["coverage_fraction"], sentinel_available=True,
    )
    active_candidate.update({"psi": psi_result["psi"], "psi_status": psi_result["status"],
                            "psi_reasons": psi_result["reasons"], "psi_weights": psi_result.get("weights")})

st.markdown('<div class="panel"><div class="panel-title">Explainable similarity and mapped composition</div><div class="panel-sub">Terrain and geology components remain visible as independent scores</div>', unsafe_allow_html=True)
if filtered_candidates:
    candidate = active_candidate
    tcol, gcol = st.columns([.8, 1.2], gap="large")
    with tcol:
        st.metric("Terrain proxy score", f"{candidate['terrain_score']:.0f}/100")
        st.caption(f"Annular relief ratio: {candidate['annular_relief_ratio']:.3f} · {active_target['body']} depth/diameter: {target_ratio:.3f}")
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
        st.metric("Earth flank-band mean slope", f"{candidate['flank_mean_slope_deg']:.1f}°")
        st.caption("Descriptive only; not scored against the target's median interior wall slope, which uses a different statistic and spatial baseline.")
        st.metric("Planetary similarity index", "Unavailable" if candidate["psi"] is None else f"{candidate['psi']:.0f}/100")
        if candidate.get("psi_weights"):
            weights_note = " Formula: " + " + ".join(f"{weight:.0%} {key}" for key, weight in candidate["psi_weights"].items()) + "."
        elif sentinel_metric is not None:
            weights_note = "."  # PSI was withheld after an incomplete Sentinel candidate overlay.
        else:
            weights_note = " Formula: 65% terrain proxy + 35% geology score normalized from /75 to /100." if candidate["psi"] is not None else "."
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
    st.markdown(f"""**Scoring recipe**  
Terrain proxy score = max(0, 100 × (1 − |Earth annular-relief ratio − target depth/diameter| / target depth/diameter)). Earth annular relief is the mean elevation in the outer 0.75–1.00 radius ring minus the mean elevation in the inner 0–0.25 radius disk. It is **not** a measured crater depth. The Earth flank-band slope is descriptive and excluded because it is not methodologically comparable to the target median wall slope.  
{geology_recipe}  
PSI first normalizes the 0–75 GLiM geology score to 0–100. Without Sentinel, PSI = 0.65 × terrain proxy + 0.35 × normalized geology. With a co-registered Sentinel-2 chip covering the selected footprint, PSI = 0.55 × terrain + 0.25 × normalized geology + 0.20 × spectral heterogeneity. PSI is withheld if GLiM coverage is below 90%, valid DEM coverage below 95%, the footprint spans fewer than 6 DEM pixels, Sentinel coverage (when used) is below 90%, or overlapping GLiM polygons exceed 1%. These are screening gates, not a calibrated uncertainty model.""")
else:
    st.info("Run the regional candidate scan, or lower the current coverage filter, to inspect terrain and GLiM scores.")
st.markdown('</div>', unsafe_allow_html=True)

with st.expander("Candidate comparison and ranking diagnostics", expanded=False):
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
        diagnostic_candidates = [item for item in ranking_pool
                                 if item["dem_coverage_fraction"] >= minimum_coverage
                                 and item["geology_coverage_fraction"] >= minimum_coverage]
        diagnostic_rows = []
        for label, wt, wg, kind in [("Terrain-only ablation", 1.0, 0.0, "terrain"),
                                    ("Geology-only ablation", 0.0, 1.0, "geology"),
                                    (f"Sensitivity mix ({terrain_weight:.0%}/{1-terrain_weight:.0%})", terrain_weight, 1-terrain_weight, "blend"),
                                    ("Baseline PSI (65/35)", 0.65, 0.35, "baseline")]:
            scored = []
            for item in diagnostic_candidates:
                geology_score = item["geology_score"]
                if kind == "terrain":
                    score = item["terrain_score"]
                elif kind == "geology":
                    score = None if geology_score is None else geology_score * (100 / 75)
                elif kind == "baseline":
                    score = item["psi"]
                else:
                    normalized_geology = None if geology_score is None else geology_score * (100 / 75)
                    score = (wt * item["terrain_score"] + wg * normalized_geology
                             if item["psi_eligible"] and geology_score is not None else None)
                scored.append((item, score))
            scored.sort(key=lambda pair: pair[1] if pair[1] is not None else -1, reverse=True)
            for rank, (item, score) in enumerate(scored, 1):
                diagnostic_rows.append({"Ranking scenario": label, "Rank": rank,
                                        "Candidate": f"{item['tile']} · {format_coordinates(item['lat'], item['lon'], 4)}",
                                        "Score": score, "PSI eligible": item["psi_eligible"]})
        st.dataframe(pd.DataFrame(diagnostic_rows), hide_index=True, use_container_width=True)
        st.caption(f"Diagnostics cover {len(diagnostic_candidates)} terrain-stratified pool samples; the map shows the top spatially separated results. Terrain-only and geology-only rows are component ablations, not PSI. Sensitivity mixes and baseline PSI are withheld wherever the base evidence gates fail. These diagnostics show ranking dependence, not independent validation or calibrated accuracy.")
    else:
        st.info("Ranking diagnostics become available after candidates pass the current coverage filter.")
    st.warning("No independent, preselected terrestrial analog validation set is bundled. This prototype reports ranking diagnostics and documented limitations; it does not claim measured validation accuracy.")

render_provenance(st)

render_mobility(st, pd, np, go, particle_size_summary, bekker_wong_scenario)
