"""Render provenance UI feature."""


def render_provenance(st):
    with st.expander("Data provenance and product metadata"):
        st.markdown("""
    | Input | Role and source status | Version / license metadata |
    | --- | --- | --- |
    | LROC NAC PDS EDR and supplied previews | NASA challenge/project observation products; DN statistics are image brightness only | Product IDs and acquisition metadata are read from PDS labels; license/version details are not present in the local manifest |
    | Copernicus GLO-30 GeoTIFFs | Supplementary Earth elevation input for the Columbia regional scan | Local filenames identify product family; specific release/license metadata are not included in the provided tiles |
    | GLiM clipped GeoPackage | Supplementary lithology polygons | Source-map scale and release/license metadata are absent from the clipped file |
    | Published McGetchin measurements / NAC DTM link | Supplementary published reference and optional lunar topography | Citation/link shown with the profile uploader; local source file attribution retained |
    | Sentinel-2 L2A | Optional user-supplied spectral input | Product identifiers, processing baseline, scale and license must be checked from the uploaded product metadata by the user |

    The scan is regional and uses local cached files; it does not download source data during the scan. The LROC previews are source imagery, while the elevation window, candidate map, slope histogram, contact mask and scoring plots are derived visualizations.
    """)
