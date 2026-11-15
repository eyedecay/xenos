"""Render candidate panel UI feature."""

def render_candidate_panel(st, right, go, np, math, filtered_candidates, candidates, format_coordinates, process_sentinel2_uploads, sentinel_candidate_metrics):
    candidate = None
    sentinel_metric = None
    with right:
        st.markdown('<div class="panel"><div class="panel-title">Earth analog viewer</div><div class="panel-sub">Copernicus GLO-30 candidate topography · optional Sentinel-2 L2A</div>', unsafe_allow_html=True)
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
            lons, lats, custom, colors, labels = [], [], [], [], []
            for view_index, item in enumerate(filtered_candidates):
                lons.append(item["lon"]); lats.append(item["lat"]); custom.append(view_index)
                score = item["psi"] if item["psi"] is not None else item["terrain_score"]
                colors.append(score)
                classes = ", ".join(row["litho"] for row in item["geology"]["composition"][:3]) or "no mapped class"
                labels.append(f"{view_index + 1}. {item['tile']} · PSI {item['psi'] if item['psi'] is not None else 'withheld'} · GLiM: {classes} · coverage DEM {item['dem_coverage_fraction']:.0%}/GLiM {item['geology']['coverage_fraction']:.0%}")
                # Approximate circular footprint in geographic degrees for regional overview only.
                radius_lat = item["diameter_m"] / 2 / 111320
                radius_lon = radius_lat / max(math.cos(math.radians(item["lat"])), 0.01)
                theta = np.linspace(0, 2 * np.pi, 49)
                map_fig.add_trace(go.Scatter(x=item["lon"] + radius_lon * np.cos(theta),
                                             y=item["lat"] + radius_lat * np.sin(theta),
                                             mode="lines", line=dict(color="#c5f36b", width=1),
                                             hoverinfo="skip", showlegend=False))
            map_fig.add_trace(go.Scatter(x=lons, y=lats, mode="markers", customdata=custom,
                                         text=labels, hovertemplate="%{text}<br>%{y:.4f}, %{x:.4f}<extra></extra>",
                                         marker=dict(size=11, color=colors, colorscale="Viridis",
                                                     cmin=0, cmax=100, showscale=True,
                                                     colorbar=dict(title="PSI or terrain", tickfont=dict(color="#9eaaa9"))),
                                         name="Candidates"))
            map_fig.update_layout(height=300, margin=dict(l=4,r=4,t=4,b=4),
                                  paper_bgcolor="#11191d", plot_bgcolor="#11191d",
                                  font_color="#9eaaa9", xaxis_title="Longitude (°)",
                                  yaxis_title="Latitude (°)", clickmode="event+select")
            map_event = st.plotly_chart(map_fig, use_container_width=True,
                                        config={"displayModeBar": False}, on_select="rerun",
                                        selection_mode="points", key="candidate_map")
            st.caption("Map shows retained candidate locations and approximate circular footprints. The dashed rectangle is the tile bounding extent, not an authoritative search-region boundary; colors show PSI when available and terrain score when PSI is withheld.")
            try:
                selected_points = map_event.selection.points
                if selected_points:
                    st.session_state["candidate_view_index"] = int(selected_points[0]["customdata"])
            except (AttributeError, KeyError, IndexError, TypeError, ValueError):
                pass
            labels = []
            for i, item in enumerate(filtered_candidates):
                psi_label = "n/a" if item["psi"] is None else f"{item['psi']:.0f}"
                labels.append(f"{i+1:02d} · {item['tile']} · base PSI {psi_label} · {format_coordinates(item['lat'], item['lon'])}")
            if st.session_state.get("candidate_view_index", 0) >= len(filtered_candidates):
                st.session_state["candidate_view_index"] = 0
            index = st.selectbox("Candidate footprint", range(len(labels)), format_func=lambda i: labels[i],
                                 key="candidate_view_index")
            candidate = filtered_candidates[index]
            fig = go.Figure(go.Heatmap(z=candidate["window"], colorscale=[[0,"#172a32"],[.35,"#315b50"],[.65,"#a9b66a"],[1,"#e4c77b"]], colorbar=dict(title="m", tickfont=dict(color="#95a6a5")), hovertemplate="row %{y}<br>col %{x}<br>%{z:.1f} m<extra></extra>"))
            fig.update_layout(height=300, margin=dict(l=0,r=0,t=4,b=0), paper_bgcolor="#11191d", plot_bgcolor="#11191d", font_color="#9eaaa9", xaxis_title="Circular footprint samples", yaxis_title="Circular footprint samples")
            st.plotly_chart(fig, use_container_width=True, config={"displayModeBar":False})
            st.caption(f"Center {format_coordinates(candidate['lat'], candidate['lon'], 5)} · circular footprint {candidate['diameter_m']} m · annular relief {candidate['annular_relief_m']:.1f} m (not crater depth)")
            st.markdown(f'<span class="tag">COPERNICUS DEM</span><span class="tag">GLiM: {len(candidate["geology"]["composition"])} mapped classes</span><span class="tag">NOT A CONFIRMED CRATER</span>', unsafe_allow_html=True)
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
            if candidates:
                st.info("No discovered candidates meet the current DEM and GLiM coverage filter. Lower the sidebar threshold to inspect them.")
            else:
                st.info("No positive-relief windows passed the available-data screen. Adjust the target values or check scan dependencies.")
        st.markdown('</div>', unsafe_allow_html=True)
    return candidate, sentinel_metric
