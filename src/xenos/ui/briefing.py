"""Readable, evidence-first summary for the selected candidate."""


def render_candidate_briefing(st, candidate):
    if candidate is None:
        return
    geology = candidate["geology"]
    top_classes = ", ".join(item["litho"] for item in geology.get("composition", [])[:3]) or "no mapped GLiM class"
    st.markdown("### Selected Earth candidate")
    st.markdown(
        f"**{candidate.get('region', 'Earth region')}** at {candidate['lat']:.4f}, {candidate['lon']:.4f}."
    )
    st.caption(
        f"Mapped classes: {top_classes}. Evidence coverage: DEM {candidate['dem_coverage_fraction']:.0%}, "
        f"GLiM {geology['coverage_fraction']:.0%}. This is a sampled Earth analog candidate; "
        "it does not establish a lunar/Martian equivalent or landing-site safety."
    )
    if candidate.get("psi_reasons"):
        st.warning("PSI withheld or qualified: " + " · ".join(candidate["psi_reasons"]))
