"""Render planetary target provenance, imagery, and measured DTM context."""


def render_reference_panel(st, left, is_default_target, REFERENCE, MORPH, LUNAR, go,
                           characterize_lroc_products, active_target, target_dtm_profile,
                           target_dtm_metrics, target_metric_name):
    with left:
        title = "Lunar target · imagery and elevation" if is_default_target else "Planetary reference · cited crater profile"
        subtitle = "McGetchin reference + LROC NAC observations" if is_default_target else "Cited target profile and optional planetary DTM"
        st.markdown(f'<div class="panel"><div class="panel-title">{title}</div><div class="panel-sub">{subtitle}</div>', unsafe_allow_html=True)
        if is_default_target:
            st.markdown(f'<span class="tag">{REFERENCE["target"]["name"]}</span><span class="tag">{REFERENCE["target"]["body"]}</span><span class="tag">PUBLISHED STUDY CROSS-CHECK</span>', unsafe_allow_html=True)
            st.caption(f"Diameter {MORPH['diameter_m']} m · published depth {MORPH['average_depth_m']} m · interior slope median / max {MORPH['average_wall_slope_deg']}° / {MORPH['maximum_wall_slope_deg']}° · rim height {MORPH['average_rim_height_m']} m")
            st.caption("The bundled JSON still says ‘To be verified’; these values and the mare–highland setting were cross-checked against the published study. The source file is unchanged.")
            image_phase = st.radio("LROC NAC context", ["Before impact", "After impact"], horizontal=True)
            preview = LUNAR / ("Pre-impact.png" if image_phase == "Before impact" else "Post-impact.png")
            st.image(str(preview), caption=f"{image_phase} · bundled LROC NAC target context", use_container_width=True)
        else:
            st.markdown(f"**{active_target['name']}** · {active_target['body']} · crater diameter {active_target['diameter_m']:g} m · published depth {active_target['depth_m']:g} m")
            st.caption(f"Source citation supplied by user: {active_target['source_citation']}.")
        if target_dtm_profile:
            profile = target_dtm_profile
            profile_fig = go.Figure(go.Scatter(x=profile["distance_m"], y=profile["elevation_m"], mode="lines", line=dict(color="#c5f36b", width=2)))
            profile_fig.update_layout(height=220, margin=dict(l=8,r=8,t=8,b=8), paper_bgcolor="#11191d", plot_bgcolor="#11191d", font_color="#9eaaa9", xaxis_title="North–south distance from crater center (m)", yaxis_title="Planetary DTM elevation (m)")
            st.plotly_chart(profile_fig, use_container_width=True, config={"displayModeBar":False})
            st.caption(f"Profile source: {profile['source']} · CRS {profile['crs']} · valid samples {profile['valid_fraction']:.0%}.")
        if target_dtm_metrics:
            m = target_dtm_metrics
            st.metric("Planetary DTM annular relief", f"{m['annular_relief_m']:.2f} m" if m["eligible"] else "Not eligible")
            st.caption(f"Same inner disk / outer annulus as the Earth DEM · ratio {m['annular_relief_ratio']:.4f} · footprint coverage {m['dem_coverage_fraction']:.1%} · {m['cells_across']:.1f} cells across · {m['source']}.")
            if m.get("analysis_pixel_size_m"):
                st.caption(f"Planetary source spacing: {m['native_pixel_size_m']:.1f} m/px · measurement grid: {m['analysis_pixel_size_m']:.1f} m/px. When the selected Earth DEM is coarser, the planetary DTM is average-resampled to that spacing. This improves scale comparability but does not remove differences between sensors or terrain products.")
            st.caption(f"Terrain comparison uses {target_metric_name}.")
        elif is_default_target:
            st.info("No suitable planetary DTM supplied. Terrain ranking uses the published depth/diameter ratio as an unlike-measure screening proxy; it is not a crater-depth match.")
            st.markdown("[Published lunar crater DTM and aligned terrain data (Zenodo)](https://doi.org/10.5281/zenodo.21481275)")
        else:
            st.caption(f"Earth annular relief is screened against the cited depth/diameter ratio ({target_metric_name}); no planetary elevation is inferred.")
        if is_default_target:
            try:
                products = characterize_lroc_products()
                with st.expander("PDS observation metadata and image statistics"):
                    for product in products:
                        st.markdown(f"**{product['phase']} — {product['product_id']}**")
                        st.write(f"{product['instrument']} · {product['shape']} pixels · producer {product['producer']} · {product['start_time']} to {product['stop_time']} UTC")
                        st.write("Centered-crop EDR brightness statistics (digital numbers; not elevation):")
                        st.json(product["stats"])
            except Exception as exc:
                st.warning(f"Could not read the PDS EDR products: {exc}")
                st.image(str(LUNAR / "Post-impact.png"), caption="Bundled post-impact preview (not used for terrain metrics)", use_container_width=True)
            st.markdown('<p class="smallnote">LROC EDR brightness values are not elevation. Only the uploaded body-referenced DTM is used to derive the planetary annular-relief comparison metric.</p></div>', unsafe_allow_html=True)
        else:
            st.markdown("</div>", unsafe_allow_html=True)
