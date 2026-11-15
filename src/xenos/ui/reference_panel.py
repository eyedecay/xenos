"""Render reference panel UI feature."""

def render_reference_panel(st, left, is_default_target, REFERENCE, MORPH, LUNAR, go, lunar_dtm_profile, characterize_lroc_products, active_target):
    with left:
        if is_default_target:
            st.markdown('<div class="panel"><div class="panel-title">Lunar target · imagery and elevation</div><div class="panel-sub">McGetchin reference values + LROC NAC observations</div>', unsafe_allow_html=True)
            st.markdown(f'<span class="tag">{REFERENCE["target"]["name"]}</span><span class="tag">{REFERENCE["target"]["body"]}</span><span class="tag">PUBLISHED STUDY CROSS-CHECK</span>', unsafe_allow_html=True)
            st.caption(f"Diameter {MORPH['diameter_m']} m · depth {MORPH['average_depth_m']} m · interior slope median / max {MORPH['average_wall_slope_deg']}° / {MORPH['maximum_wall_slope_deg']}° · rim height {MORPH['average_rim_height_m']} m")
            st.caption("The bundled JSON still says ‘To be verified’; these values and the mare–highland setting were cross-checked against the published study. The source file is unchanged.")
            image_phase = st.radio("LROC NAC context", ["Before impact", "After impact"], horizontal=True)
            preview = LUNAR / ("Pre-impact.png" if image_phase == "Before impact" else "Post-impact.png")
            st.image(str(preview), caption=f"{image_phase} · bundled LROC NAC target preview", use_container_width=True)
            st.markdown("**Lunar elevation profile**")
            lunar_dtm_upload = st.file_uploader("Upload the published NAC lunar DTM GeoTIFF", type=["tif", "tiff"], key="lunar_dtm")
            if lunar_dtm_upload:
                try:
                    profile = lunar_dtm_profile(lunar_dtm_upload.getvalue())
                    profile_fig = go.Figure(go.Scatter(x=profile["distance_m"], y=profile["elevation_m"], mode="lines", line=dict(color="#c5f36b", width=2)))
                    profile_fig.update_layout(height=220, margin=dict(l=8,r=8,t=8,b=8), paper_bgcolor="#11191d", plot_bgcolor="#11191d", font_color="#9eaaa9", xaxis_title="North–south distance from crater center (m)", yaxis_title="Lunar DTM elevation (m)")
                    st.plotly_chart(profile_fig, use_container_width=True, config={"displayModeBar":False})
                    st.caption(f"Profile source: {profile['source']} · CRS {profile['crs']} · valid samples {profile['valid_fraction']:.0%}.")
                except Exception as exc:
                    st.warning(f"Could not profile this lunar DTM: {exc}")
            else:
                st.info("No lunar DTM is bundled. Upload the NAC stereo DTM to plot a measured profile; this is a better-matched crater-scale product than LOLA grids.")
                st.markdown("[Published crater DTM and aligned terrain data (Zenodo)](https://doi.org/10.5281/zenodo.21481275)")
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
            st.markdown('<p class="smallnote">The supplied before/after previews show the lunar target. PDS DN statistics describe image values only. Elevation profiles require the published lunar DTM upload; no profile is inferred from NAC brightness or from a synthetic curve.</p></div>', unsafe_allow_html=True)
        else:
            st.markdown('<div class="panel"><div class="panel-title">Planetary reference · cited crater profile</div><div class="panel-sub">No mission imagery or elevation data bundled for this uploaded target</div>', unsafe_allow_html=True)
            st.markdown(f"**{active_target['name']}** · {active_target['body']} · crater diameter {active_target['diameter_m']:g} m · depth {active_target['depth_m']:g} m")
            st.caption(f"Source citation supplied by user: {active_target['source_citation']}. Earth screening uses the cited depth/diameter ratio; no planetary imagery is inferred.")
            st.markdown('</div>', unsafe_allow_html=True)
