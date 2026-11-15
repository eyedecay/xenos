"""Optional independent analog-site validation controls."""
from ..validation import validate_known_sites
import pandas as pd


def render_validation(st, ranking_pool):
    st.markdown("**Optional known-site validation**")
    upload = st.file_uploader(
        "Upload independent analog-site coordinates CSV",
        type=["csv"], key="validation_sites",
        help="Required columns: site_name, latitude, longitude. Optional tolerance_m. Keep these sites out of score/profile tuning.",
    )
    if upload is None:
        st.caption("No reference validation set is bundled. The app makes no accuracy claim until independent sites are provided.")
        st.download_button("Download validation CSV template",
                           "site_name,latitude,longitude,tolerance_m\n",
                           "xenos_validation_template.csv", "text/csv")
        return
    radius = st.number_input("Default spatial match tolerance (m)", min_value=100.0,
                             max_value=100000.0, value=5000.0, step=500.0,
                             key="validation_tolerance")
    try:
        rows, summary = validate_known_sites(upload.getvalue(), ranking_pool, radius)
        st.metric("Known sites within tolerance", f"{summary['sites_with_candidate_within_tolerance']}/{summary['sites']}")
        st.caption(f"Median nearest sampled-candidate distance: {summary['median_nearest_distance_m']:.0f} m. These diagnostics depend on candidate sampling density and the supplied tolerance; they are not calibrated accuracy.")
        result_frame = pd.DataFrame(rows)
        st.dataframe(result_frame, hide_index=True, use_container_width=True)
        st.download_button("Download validation results CSV", result_frame.to_csv(index=False),
                           "xenos_validation_results.csv", "text/csv")
    except Exception as exc:
        st.error(f"Validation could not be calculated: {exc}")
