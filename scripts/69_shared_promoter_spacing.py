#!/usr/bin/env python
"""69_shared_promoter_spacing.py — how close are the two TSSs of a pair?

A shared bidirectional promoter places two initiation sites a short, fixed
distance apart, typically on the order of a nucleosome-free region, rather than
anywhere within the 1 kb that defines the class. FANTOM5 calls peaks per strand,
so two genes never share a peak identifier even when they share a promoter. The
informative quantity is therefore the spacing itself.

The null is the spacing between protein-coding TSSs that face each other on
opposite strands, which is what the genome produces without any lncRNA. If
divergent pairs sit at the same spacing as PCG-PCG pairs, the lncRNA is using
the ordinary bidirectional arrangement. If they sit further apart, the class
window rather than promoter sharing is doing the work.

Outputs
  results/partner_paired/tss_spacing.tsv
  results/partner_paired/tss_spacing_pairs.tsv
  results/step_summaries/69_shared_promoter_spacing.md
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "partner_paired")
WIN = 1000


def nearest_opposite(tss_df, exclude_same_gene=True):
    """Distance from every TSS to the nearest opposite-strand TSS within WIN."""
    out = []
    for chrom, g in tss_df.groupby("chrom"):
        plus = np.sort(g.tss[g.strand == "+"].values)
        minus = np.sort(g.tss[g.strand == "-"].values)
        if not len(plus) or not len(minus):
            continue
        for s, arr, other in (("+", plus, minus), ("-", minus, plus)):
            j = np.searchsorted(other, arr)
            for k, pos in enumerate(arr):
                best = np.inf
                for jj in (j[k] - 1, j[k], j[k] + 1):
                    if 0 <= jj < len(other):
                        best = min(best, abs(other[jj] - pos))
                if np.isfinite(best) and best <= WIN:
                    out.append(best)
    return np.array(out, dtype=float)


def main():
    with step("69_shared_promoter_spacing") as st:
        c = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                        low_memory=False)
        c["tss"] = pd.to_numeric(c.tss, errors="coerce")
        c = c.dropna(subset=["tss"])

        div = c[(c.cls == "divergent") & c.partner_pcg.notna()].copy()
        partner_tss = c.set_index("gene_id").tss
        div["partner_tss"] = div.partner_pcg.map(partner_tss)
        div = div.dropna(subset=["partner_tss"])
        div["spacing"] = (div.tss - div.partner_tss).abs()
        st.note("divergent pairs={} median spacing={} bp (IQR {}–{})".format(
            len(div), int(div.spacing.median()),
            int(div.spacing.quantile(.25)), int(div.spacing.quantile(.75))))

        pcg = c[c.group == "pcg"][["chrom", "tss", "strand"]]
        null = nearest_opposite(pcg)
        st.note("PCG-PCG opposite-strand spacing: n={} median={} bp "
                "(IQR {}–{})".format(len(null), int(np.median(null)),
                                    int(np.percentile(null, 25)),
                                    int(np.percentile(null, 75))))

        w, p = mannwhitneyu(div.spacing, null, alternative="two-sided")
        rows = []
        for q in (0.1, 0.25, 0.5, 0.75, 0.9):
            rows.append(dict(quantile=q,
                             divergent=float(div.spacing.quantile(q)),
                             pcg_null=float(np.quantile(null, q))))
        T = pd.DataFrame(rows)
        T.to_csv(safe_write_path(os.path.join(OUT, "tss_spacing.tsv")),
                 sep="\t", index=False)
        div[["gene_id", "partner_pcg", "spacing", "div_dist_opp"]].to_csv(
            safe_write_path(os.path.join(OUT, "tss_spacing_pairs.tsv")),
            sep="\t", index=False)

        frac_short = float((div.spacing <= 100).mean())
        frac_null_short = float((null <= 100).mean())
        st.note("spacing <= 100 bp: {:.1%} of pairs versus {:.1%} of the PCG "
                "null; Mann-Whitney P={:.3g}".format(
                    frac_short, frac_null_short, p))
        T2 = pd.DataFrame([dict(
            n_pairs=len(div), n_null=len(null),
            median_divergent=float(div.spacing.median()),
            median_pcg_null=float(np.median(null)),
            frac_le_100_divergent=frac_short,
            frac_le_100_null=frac_null_short,
            frac_le_250_divergent=float((div.spacing <= 250).mean()),
            frac_le_250_null=float((null <= 250).mean()),
            mannwhitney_p=float(p))])
        T2.to_csv(safe_write_path(os.path.join(
            OUT, "tss_spacing_summary.tsv")), sep="\t", index=False)
        st.summary_md(os.path.join(RES, "step_summaries",
                                   "69_shared_promoter_spacing.md"),
                      "69 spacing between the two TSSs of a divergent pair",
                      df_to_md(T, ".1f") + "\n\n" + df_to_md(T2, ".4f"))


if __name__ == "__main__":
    main()
