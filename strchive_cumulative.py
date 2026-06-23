#!/usr/bin/env python3
"""Cumulative number of known repeat-expansion disorder genes over time,
built from STRchive-loci.json (gene + year), filtered by evidence level.

Reproduces the style of Depienne et al. "Growing number of repeat expansion
disorders" figure.
"""
import json
import re
import argparse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

KEEP_EVIDENCE = {"Definitive", "Strong", "Moderate", "Limited", "Provisional"}


# Manual overrides for loci with missing/ambiguous location_in_gene
CLASS_OVERRIDES = {
    "GOLGA8A": "intron",  # location_in_gene is null in STRchive; expansion is intronic
}


def classify(loc, gene=None):
    """Region class from the location_in_gene string."""
    if gene in CLASS_OVERRIDES:
        return CLASS_OVERRIDES[gene]
    if not loc:
        return "other"
    low = loc.lower()
    if "promoter" in low:
        return "promoter"
    if "utr" in low:
        return "UTR"
    if "intron" in low:
        return "intron"
    if "lncrna" in low or "non-coding" in low or "noncoding" in low:
        return "other"
    if "coding" in low or "exon" in low:
        return "coding"
    return "other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default="STRchive-loci.json")
    ap.add_argument("--out", default="strchive_cumulative.png")
    args = ap.parse_args()

    data = json.load(open(args.json))
    loci = [
        x for x in data
        if set(x.get("evidence") or []) & KEEP_EVIDENCE and x.get("year")
    ]
    for x in loci:
        x["year"] = int(re.search(r"\d{4}", str(x["year"])).group())
    # one point per locus, ordered by year; break ties by gene name
    loci.sort(key=lambda x: (x["year"], x["gene"]))

    colors = {
        "coding": "#1f77b4",
        "UTR": "#ff7f0e",
        "promoter": "#2ca02c",
        "intron": "#d62728",
        "other": "#7f7f7f",
    }
    markers = {
        "coding": "s",
        "UTR": "o",
        "promoter": "D",
        "intron": "^",
        "other": "v",
    }

    fig, ax = plt.subplots(figsize=(14, 9))

    # long-read era shading (open-ended: extends to the right plot edge)
    ax.axvspan(2018.5, max(x["year"] for x in loci) + 100, color="#cfe2f3", alpha=0.5, zorder=0)

    # No month is available in STRchive, so loci sharing a year are spread
    # evenly across that year's width (cosmetic only) to avoid vertical stacks.
    from collections import Counter as _C
    year_counts = _C(x["year"] for x in loci)
    seen = {}

    xs, ys = [], []
    for i, x in enumerate(loci, start=1):
        cls = classify(x.get("location_in_gene"), x["gene"])
        n = year_counts[x["year"]]
        k = seen.get(x["year"], 0)
        seen[x["year"]] = k + 1
        xpos = x["year"] + (k + 0.5) / n * 0.85 if n > 1 else x["year"]
        ax.scatter(xpos, i, s=45, color=colors[cls], marker=markers[cls],
                   zorder=3, edgecolors="none")
        label = x["gene"]
        ax.annotate(label, (xpos, i), xytext=(4, -2),
                    textcoords="offset points", fontsize=8, va="center")
        xs.append(xpos)
        ys.append(i)

    ax.plot(xs, ys, color="grey", lw=0.6, zorder=2, alpha=0.6)

    ax.set_xlabel("Year", fontsize=15)
    ax.set_ylabel("Cumulative number of repeat-expansion loci / disorders", fontsize=15)
    ax.tick_params(axis="both", labelsize=12)
    ax.set_title("Growing number of repeat expansion disorders", fontsize=18, fontweight="bold")
    ax.text(2022.5, 8, "Long-read\nsequencing\nera", color="firebrick",
            fontsize=13, ha="center", fontweight="bold")

    legend = [
        Line2D([0], [0], marker=markers[c], color="w", markerfacecolor=colors[c],
               markersize=9, label=c)
        for c in ["coding", "UTR", "promoter", "intron", "other"]
    ]
    ax.legend(handles=legend, loc="upper left", frameon=True, fontsize=12)

    ax.set_xlim(min(xs) - 1, max(xs) + 2)
    ax.set_ylim(0, len(loci) + 2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    from collections import Counter
    counts = Counter(classify(x.get("location_in_gene"), x["gene"]) for x in loci)
    print(f"Wrote {args.out} with {len(loci)} loci: " +
          ", ".join(f"{k}={v}" for k, v in counts.most_common()))


if __name__ == "__main__":
    main()
