"""Compare ranking scenarios on the same PSI-eligible sample."""
def compare_rankings(pool, terrain_weight, k=10):
    eligible = [item for item in pool if item.get("psi") is not None]
    baseline = sorted(range(len(eligible)), key=lambda i: eligible[i]["psi"], reverse=True)
    shortlist = set(baseline[:k])
    summaries, rankings = [], []
    for label, weight in [("Baseline PSI (65/35)", None), ("Terrain only", 1.0),
                          ("Geology only", 0.0), (f"Selected mix ({terrain_weight:.0%} terrain)", terrain_weight)]:
        scores = [c["psi"] if weight is None else weight * c["terrain_score"] +
                  (1 - weight) * (c["geology_score"] * (100 / 75)) for c in eligible]
        order = sorted(range(len(scores)), key=scores.__getitem__, reverse=True)
        shared = len(shortlist.intersection(order[:k]))
        tied = len(order) > k and scores[order[k - 1]] == scores[order[k]]
        summaries.append({"Scenario": label, "Samples compared": len(order), "Shortlist size": len(shortlist),
                          "Shared with baseline": shared, "Tie at cutoff": tied})
        rankings.extend((label, rank, eligible[i], scores[i]) for rank, i in enumerate(order, 1))
    return summaries, rankings
