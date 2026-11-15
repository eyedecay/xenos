"""Download buttons for reproducible candidate results."""
from ..outputs import candidates_csv, candidates_geojson, run_metadata_json


def render_exports(st, candidates, metadata):
    with st.expander("Export results and reproduce this scan", expanded=False):
        if not candidates:
            st.info("Run a scan to export candidate results.")
            return
        st.caption("Exports include retained candidate metrics and a run manifest with target citation, input fingerprints, score formula, thresholds, and scan summary. Large source rasters and the candidate elevation arrays are not embedded.")
        columns = st.columns(3)
        columns[0].download_button("Download candidate CSV", candidates_csv(candidates),
                                  "xenos_candidates.csv", "text/csv", use_container_width=True)
        columns[1].download_button("Download GeoJSON", candidates_geojson(candidates),
                                  "xenos_candidates.geojson", "application/geo+json", use_container_width=True)
        columns[2].download_button("Download run manifest", run_metadata_json(metadata),
                                  "xenos_run_manifest.json", "application/json", use_container_width=True)
