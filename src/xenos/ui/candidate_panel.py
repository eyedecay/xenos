"""Render candidate panel UI feature."""

def render_candidate_panel(st, right, go, np, math, filtered_candidates, candidates, format_coordinates, process_sentinel2_uploads, sentinel_candidate_metrics):
    candidate = None
    sentinel_metric = None
    with right:
        st.markdown('<div class="panel"><div class="panel-title">Earth candidate map</div><div class="panel-sub">Circular candidate footprints · Copernicus GLO-30 terrain</div>', unsafe_allow_html=True)
        if st.session_state.scan_error:
            st.error(f"Candidate scan could not complete: {st.session_state.scan_error}")
            st.caption("Install the dependencies in requirements.txt, then restart the app.")
        elif filtered_candidates:
            map_fig = go.Figure()
            if st.session_state.scan_summary and st.session_state.scan_summary.get("region_extent_lonlat"):
                min_lon, min_lat, max_lon, max_lat = st.session_state.scan_summary["region_extent_lonlat"]
                map_fig.add_shape(type="rect", x0=min_lon, y0=min_lat, x1=max_lon, y1=max_lat,
                                  line=dict(color="#71837f", dash="dash", width=1),
                                  fillcolor="rgba(0,0,0,0)")
            eligible_points, withheld_points = [], []
            for view_index, item in enumerate(filtered_candidates):
                score = item["psi"]
                classes = ", ".join(row["litho"] for row in item["geology"]["composition"][:3]) or "no mapped class"
                score_label = f"PSI {score:.0f}/100" if score is not None else "PSI withheld"
                label = (f"{view_index + 1}. {item['tile']} · {score_label} · GLiM: {classes} · "
                         f"coverage DEM {item['dem_coverage_fraction']:.0%}/GLiM {item['geology']['coverage_fraction']:.0%}")
                point = (item["lon"], item["lat"], view_index, label)
                (eligible_points if score is not None else withheld_points).append((point, score))
                # Approximate circular footprint in geographic degrees for regional overview only.
                radius_lat = item["diameter_m"] / 2 / 111320
                radius_lon = radius_lat / max(math.cos(math.radians(item["lat"])), 0.01)
                theta = np.linspace(0, 2 * np.pi, 49)
                map_fig.add_trace(go.Scatter(x=item["lon"] + radius_lon * np.cos(theta),
                                             y=item["lat"] + radius_lat * np.sin(theta),
                                             mode="lines", line=dict(color="#c5f36b", width=1),
                                             hoverinfo="skip", showlegend=False))
            if eligible_points:
                map_fig.add_trace(go.Scatter(
                    x=[p[0][0] for p in eligible_points], y=[p[0][1] for p in eligible_points],
                    mode="markers", customdata=[p[0][2] for p in eligible_points],
                    text=[p[0][3] for p in eligible_points],
                    hovertemplate="%{text}<br>%{y:.4f}, %{x:.4f}<extra></extra>",
                    marker=dict(size=11, color=[p[1] for p in eligible_points], colorscale="Viridis",
                                                     cmin=0, cmax=100, showscale=True,
                                colorbar=dict(title="PSI", tickfont=dict(color="#9eaaa9"))),
                    name="PSI eligible"))
            if withheld_points:
                map_fig.add_trace(go.Scatter(
                    x=[p[0][0] for p in withheld_points], y=[p[0][1] for p in withheld_points],
                    mode="markers", customdata=[p[0][2] for p in withheld_points],
                    text=[p[0][3] for p in withheld_points],
                    hovertemplate="%{text}<br>%{y:.4f}, %{x:.4f}<extra></extra>",
                    marker=dict(size=11, color="#788486", symbol="circle-open"),
                    name="PSI withheld"))
            map_fig.update_layout(height=300, margin=dict(l=4,r=4,t=4,b=4),
                                  paper_bgcolor="#11191d", plot_bgcolor="#11191d",
                                  font_color="#9eaaa9", xaxis_title="Longitude (°)",
                                  yaxis_title="Latitude (°)")
            st.plotly_chart(map_fig, use_container_width=True,
                            config={"displayModeBar": False}, key="candidate_map")
            st.caption("Hover over markers for candidate details. Use the candidate selector below to inspect a footprint. Colors show PSI; open gray markers indicate withheld PSI. The dashed rectangle is a tile extent, not an authoritative search boundary.")
            labels = []
            for i, item in enumerate(filtered_candidates):
                psi_label = "withheld" if item["psi"] is None else f"{item['psi']:.0f}/100"
                labels.append(f"{i+1:02d} · {item['tile']} · PSI {psi_label} · {format_coordinates(item['lat'], item['lon'])}")
            if st.session_state.get("candidate_selector", 0) >= len(filtered_candidates):
                st.session_state["candidate_selector"] = 0
            index = st.selectbox("Candidate footprint", range(len(labels)), format_func=lambda i: labels[i],
                                 key="candidate_selector")
            candidate = filtered_candidates[index]
            fig = go.Figure(go.Heatmap(z=candidate["window"], colorscale=[[0,"#172a32"],[.35,"#315b50"],[.65,"#a9b66a"],[1,"#e4c77b"]], colorbar=dict(title="m", tickfont=dict(color="#95a6a5")), hovertemplate="row %{y}<br>col %{x}<br>%{z:.1f} m<extra></extra>"))
            fig.update_layout(height=300, margin=dict(l=0,r=0,t=4,b=0), paper_bgcolor="#11191d", plot_bgcolor="#11191d", font_color="#9eaaa9", xaxis_title="Circular footprint samples", yaxis_title="Circular footprint samples")
            st.plotly_chart(fig, use_container_width=True, config={"displayModeBar":False})
            st.caption(f"Center {format_coordinates(candidate['lat'], candidate['lon'], 5)} · circular footprint {candidate['diameter_m']} m · annular relief {candidate['annular_relief_m']:.1f} m (not crater depth)")
            st.markdown(f'<span class="tag">COPERNICUS DEM</span><span class="tag">GLiM: {len(candidate["geology"]["composition"])} mapped classes</span><span class="tag">NOT A CONFIRMED CRATER</span>', unsafe_allow_html=True)
            with st.expander("Optional spectral check (experimental)", expanded=False):
                st.markdown("**Sentinel-2 contact-line scan**")
                sentinel_uploads = st.file_uploader(
                    "Upload co-registered L2A reflectance chips named with B02, B03, B04, B11 and B12",
                    type=["tif", "tiff"], accept_multiple_files=True, key="sentinel_bands",
                    help="Use a small projected reflectance chip (maximum 2048 pixels/side) covering this candidate. Optional SCL band masks clouds and shadows.",
                )
                scale_col, threshold_col = st.columns(2)
                with scale_col:
                    reflectance_scale = st.number_input("Reflectance scale", min_value=0.000001, max_value=1.0,
                                                        value=1.0, format="%.6f", key="s2_scale")
                with threshold_col:
                    contact_threshold = st.number_input("Contact gradient threshold / pixel", min_value=0.001,
                                                         max_value=1.0, value=0.02, step=0.005,
                                                         key="s2_threshold")
                st.caption("Scale 1.0 is for reflectance already in 0–1 units. Raw L2A DN requires the product quantification and per-band offsets; this input applies a common multiplier only.")
                if sentinel_uploads:
                    try:
                        sentinel_data = process_sentinel2_uploads(sentinel_uploads, reflectance_scale)
                        sentinel_metric = sentinel_candidate_metrics(
                            sentinel_data, candidate["lon"], candidate["lat"], candidate["diameter_m"],
                            contact_threshold,
                        )
                        spectral_cols = st.columns(2)
                        with spectral_cols[0]:
                            st.image(sentinel_data["true_color"], caption="Sentinel-2 true color (B04/B03/B02)", use_container_width=True)
                        with spectral_cols[1]:
                            if sentinel_metric["contact_mask"] is not None:
                                st.image(sentinel_metric["contact_mask"].astype("uint8") * 255,
                                         caption="Spectral contact-gradient pixels", use_container_width=True)
                        st.metric("Sentinel spectral heterogeneity", "Unavailable" if sentinel_metric["score"] is None else f"{sentinel_metric['score']:.1f}/100")
                        st.caption(f"Candidate chip coverage {sentinel_metric['coverage_fraction']:.1%} · projected CRS {sentinel_data['crs']} · SCL {sentinel_data['scl']}. Iron proxy: normalized B04−B02; SWIR contrast: normalized B11−B12. These are screening indices, not mineral identifications.")
                        st.caption("This chip updates PSI for the selected candidate only; the regional candidate list remains ordered by DEM + GLiM until spectral coverage is available for every candidate.")
                    except Exception as exc:
                        st.error(f"Sentinel-2 processing could not complete: {exc}")
        else:
            if st.session_state.get("scan_args") is None:
                st.info("Set the reference, region, and footprint size, then run the scan to see candidates here.")
            else:
                st.info("The last scan found no positive-relief windows that passed the available-data screen. Check the selected target or region data.")
        st.markdown('</div>', unsafe_allow_html=True)
    return candidate, sentinel_metric
