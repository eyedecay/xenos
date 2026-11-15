"""Published crater and uploaded reference-site retrieval controls."""
import pandas as pd

from ..config import ROOT
from ..validation import validate_known_sites


def render_validation(st, ranking_pool, odessa_scan=False):
    st.markdown("**Check the ranking against reference sites**")
    upload = st.file_uploader(
        "Reference site CSV", type=["csv"], key="validation_sites",
        help="Supply site_name, latitude, longitude and optional tolerance_m. Choose sites and tolerances before viewing results.",
    )
    if upload is not None:
        sites = upload.getvalue()
    elif odessa_scan:
        sites = (ROOT / "examples/odessa.csv").read_bytes()
        st.caption("Odessa main crater · published location · 222 m search tolerance. This single case checks crater proximity. "
                   "[Meteoritical Bulletin](https://www.lpi.usra.edu/meteor/metbull.php?code=35275)")
    else:
        st.caption("Scan the Odessa crater study area alone for the published example, or upload sites within your scanned regions.")
        st.download_button("Download validation CSV template", "site_name,latitude,longitude,tolerance_m\n",
                           "xenos_validation_template.csv", "text/csv")
        return
    radius = (st.number_input("Default match tolerance (m)", 100.0, 100000.0, 5000.0, 500.0,
                              key="validation_tolerance") if upload is not None else 222.0)
    try:
        rows, summary = validate_known_sites(sites, ranking_pool, radius)
    except ValueError as exc:
        st.info(str(exc))
        return
    st.caption(f"{len(rows)} reference site(s) compared with {summary['eligible_pool_size']} PSI-eligible samples before map thinning. "
               "Each hit means at least one ranked sample falls within that site's tolerance.")
    comparison = [{"First K samples": min(k, summary["eligible_pool_size"]),
                   "Observed site hits": sum(row[f"hit_at_{k}"] for row in rows),
                   "Sites retrieved (%)": round(100 * summary[f"hit_at_{k}"], 1),
                   "Random expectation (%)": round(100 * summary[f"random_baseline_at_{k}"], 1)}
                  for k in (1, 5, 10)]
    st.dataframe(pd.DataFrame(comparison).drop_duplicates(), hide_index=True, use_container_width=True)
    st.caption("The baseline draws K samples uniformly without replacement from the same pool. Sampling density and match tolerance affect both results.")
    frame = pd.DataFrame(rows)
    st.dataframe(frame, hide_index=True, use_container_width=True)
    st.download_button("Download site results", frame.to_csv(index=False), "xenos_validation_results.csv", "text/csv")
