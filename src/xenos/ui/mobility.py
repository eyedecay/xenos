"""Optional PSD and terramechanics scenario UI."""

def render_mobility(st, pd, np, go, particle_size_summary, bekker_wong_scenario):
    st.markdown('<div class="panel"><div class="panel-title">Terramechanics scenario</div><div class="panel-sub">Bekker pressure-sinkage + simplified Wong-Janosi-Hanamoto traction</div>', unsafe_allow_html=True)
    st.info("No particle-size distribution or calibrated soil parameters are bundled. Upload a PSD and enter measured or deliberately chosen scenario values below; the dashboard never supplies defaults for soil properties.")
    psd_upload = st.file_uploader("Upload PSD CSV (particle_size_mm,mass_fraction)", type=["csv"], key="psd_csv")
    if psd_upload:
        try:
            psd_result = particle_size_summary(pd.read_csv(psd_upload))
            st.caption(f"PSD records: {psd_result['particle_count']} · D10 {psd_result['quantiles']['D10_mm']:.3g} mm · D50 {psd_result['quantiles']['D50_mm']:.3g} mm · D90 {psd_result['quantiles']['D90_mm']:.3g} mm. Particle sizes do not determine Bekker moduli without calibration.")
            psd_fig = go.Figure(go.Bar(x=psd_result["sizes_mm"], y=psd_result["mass_fractions"], marker_color="#6e9b83"))
            psd_fig.update_layout(height=180, margin=dict(l=8,r=8,t=4,b=8), paper_bgcolor="#11191d", plot_bgcolor="#11191d", font_color="#9eaaa9", xaxis_title="Particle size (mm)", yaxis_title="Mass fraction")
            st.plotly_chart(psd_fig, use_container_width=True, config={"displayModeBar":False})
        except Exception as exc:
            st.error(f"PSD could not be read: {exc}")

    with st.expander("Enter Bekker-Wong scenario inputs"):
        st.caption("Use one consistent SI unit system. k_c and k_phi units depend on n in p = (k_c / b + k_phi) z^n. The traction and risk outputs below are simplified screening calculations.")
        input_cols = st.columns(4)
        labels_keys = [
            ("Bekker k_c", "tm_kc"), ("Bekker k_phi", "tm_kphi"),
            ("Bekker exponent n", "tm_n"), ("Wheel width b (m)", "tm_width"),
            ("Contact length (m)", "tm_length"), ("Normal load / wheel (N)", "tm_load"),
            ("Wheel radius (m)", "tm_radius"), ("Soil cohesion c (Pa)", "tm_cohesion"),
            ("Friction angle φ (°)", "tm_phi"), ("Shear modulus K (m)", "tm_shear_k"),
            ("Shear displacement j (m)", "tm_shear_j"),
        ]
        model_values = {}
        for index_input, (label, key) in enumerate(labels_keys):
            with input_cols[index_input % len(input_cols)]:
                raw_value = st.text_input(label, key=key, placeholder="Enter measured or scenario value")
                if raw_value.strip():
                    try:
                        model_values[key] = float(raw_value)
                    except ValueError:
                        model_values[key] = None
        required_keys = [key for _, key in labels_keys]
        if any(key in model_values and model_values[key] is None for key in required_keys):
            st.error("Every populated model input must be numeric.")
        elif all(key in model_values for key in required_keys):
            try:
                mobility = bekker_wong_scenario(
                    model_values["tm_kc"], model_values["tm_kphi"], model_values["tm_n"],
                    model_values["tm_width"], model_values["tm_length"], model_values["tm_load"],
                    model_values["tm_radius"], model_values["tm_cohesion"], model_values["tm_phi"],
                    model_values["tm_shear_k"], model_values["tm_shear_j"],
                )
                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Wheel sinkage", f"{mobility['sinkage_m']*1000:.1f} mm")
                m2.metric("Drawbar pull", f"{mobility['drawbar_pull_n']:.1f} N")
                m3.metric("Drawbar efficiency", f"{mobility['drawbar_pull_efficiency']:.0%}")
                m4.metric("Mobility risk index", f"{mobility['rover_mobility_risk_index']:.0f}/100")
                curve = go.Figure(go.Scatter(x=np.array(mobility["sinkage_curve_m"])*1000,
                                             y=np.array(mobility["pressure_pa"])/1000,
                                             mode="lines", line=dict(color="#c5f36b", width=3)))
                curve.update_layout(height=230, margin=dict(l=8,r=8,t=8,b=8), paper_bgcolor="#11191d", plot_bgcolor="#11191d", font_color="#9eaaa9", xaxis_title="Sinkage (mm)", yaxis_title="Normal pressure (kPa)")
                st.plotly_chart(curve, use_container_width=True, config={"displayModeBar":False})
                st.caption(f"Pressure-sinkage uses the entered Bekker coefficients. Traction is a uniform-contact shear approximation; risk formula: {mobility['risk_method']}.")
            except Exception as exc:
                st.error(f"Scenario could not be calculated: {exc}")
        else:
            st.caption("Enter all 11 model inputs to calculate the scenario; no unmeasured values are filled automatically.")
    st.markdown('</div>', unsafe_allow_html=True)

    if st.session_state.scan_summary:
        summary = st.session_state.scan_summary
        extent = summary.get("region_extent_lonlat")
        extent_text = (f"Extent {abs(extent[0]):.3f}–{abs(extent[2]):.3f}°W, {extent[1]:.3f}–{extent[3]:.3f}°N" if extent and extent[2] < 0 else (f"Extent {extent[0]:.3f}–{extent[2]:.3f}°E, {extent[1]:.3f}–{extent[3]:.3f}°N" if extent else "extent unavailable"))
        st.caption(f"{summary['region_name']} regional scan · {extent_text} · processed {summary['tiles_read']} DEM tiles and {summary['glim_features']:,} geology polygons · screened {summary['candidate_windows']:,} windows · geology-scored {summary['candidate_pool_size']} stratified samples · retained {summary['candidate_locations']} spatially separated map candidates · PSI gates: geology ≥{summary['min_geology_coverage_for_psi']:.0%}, DEM ≥{summary['min_dem_coverage_for_psi']:.0%}, at least {summary['min_cells_across_for_psi']:.0f} pixels across.")
