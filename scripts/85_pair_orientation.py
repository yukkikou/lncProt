#!/usr/bin/env python
"""85_pair_orientation.py — head-to-head versus tail-to-tail divergent pairs.

The divergent class is defined by opposite-strand TSS proximity and does not
require the two genes to point at each other. This module splits the 2,692
pairs by gene-body geometry and reports CAGE asymmetry, window overlap and
partner-promoter overlap for each orientation, so the borrowing claim can be
read separately for the strictly head-to-head subset.

Outputs
  results/revision_checks/I_pair_orientation.tsv
  results/revision_checks/I_orientation_bins.tsv
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "revision_checks")


def main():
    with step("85_pair_orientation") as st:
        os.makedirs(OUT, exist_ok=True)
        dp = pd.read_csv(os.path.join(RES, "divergent_pairs.tsv"),
                         sep="\t").dropna()
        genes = pd.read_csv(os.path.join(WORKDIR, "data", "genes",
                                         "genes.tsv"), sep="\t").set_index(
                                             "gene_id")
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        meta["tss"] = pd.to_numeric(meta.tss, errors="coerce")
        cage = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                           low_memory=False)
        cage["base"] = cage.gene_id.astype(str).str.split(".").str[0]
        cage = cage.drop_duplicates("base").set_index("base")

        d = dp.copy()
        for c, src in (("ls", genes.strand), ("lstart", genes.start),
                       ("lend", genes.end), ("ps", genes.strand),
                       ("pstart", genes.start), ("pend", genes.end)):
            key, col = c, src
            d[key] = d.lnc_gid.map(col) if key.startswith("l") else \
                d.pcg_gid.map(col)
        d["ltss"] = d.lnc_gid.map(meta.tss)
        d["ptss"] = d.pcg_gid.map(meta.tss)
        d = d.dropna(subset=["ltss", "ptss", "lstart", "pstart"])

        # head-to-head: the two TSSs face away from each other. For a plus
        # strand partner the lncRNA TSS must lie upstream, and for a minus
        # strand partner downstream. A body-containment rule conflates this
        # with nested genes and mislabels convergent pairs.
        d["h2h"] = np.where(d.ps == "+", d.ltss < d.ptss,
                            d.ltss > d.ptss)
        d["overlap_bp"] = np.clip(
            np.minimum(d.lend, d.pend) - np.maximum(d.lstart, d.pstart),
            0, None)
        d["prom_overlap"] = np.clip(
            np.minimum(d.ltss + 200, d.ptss + 200)
            - np.maximum(d.ltss - 500, d.ptss - 500), 0, None)
        # Jaccard of the two 700 bp windows, so it is bounded by 1
        _j = d.overlap_bp / (700 * 2 - d.overlap_bp).clip(lower=1)
        d["jaccard"] = _j.clip(upper=1.0)
        d["lact"] = d.lnc_gid.astype(str).str.split(".").str[0].map(
            cage.active)
        d["pact"] = d.pcg_gid.astype(str).str.split(".").str[0].map(
            cage.active)

        rows = []
        for lab, sub in (("head_to_head", d[d.h2h]),
                         ("tail_to_tail", d[~d.h2h]),
                         ("all", d)):
            s = sub.dropna(subset=["lact", "pact"])
            a = int(((s.lact == 1) & (s.pact == 0)).sum())
            b = int(((s.lact == 0) & (s.pact == 1)).sum())
            rows.append(dict(subset=lab, n=len(sub),
                             median_spacing=float(sub.distance.median()),
                             median_body_overlap=float(
                                 sub.overlap_bp.median()),
                             median_promoter_overlap=float(
                                 sub.prom_overlap.median()),
                             median_jaccard=float(sub.jaccard.median()),
                             lnc_cage_active=float(s.lact.mean()),
                             partner_cage_active=float(s.pact.mean()),
                             lnc_only=a, partner_only=b,
                             binom_p=float(binomtest(min(a, b), a + b,
                                                     0.5).pvalue)))
        T = pd.DataFrame(rows)
        bins = pd.cut(d.distance, [0, 150, 250, 400, 600, 1000])
        R = []
        for b, sub in d.groupby(bins, observed=True):
            h = sub[sub.h2h].lact.dropna()
            t = sub[~sub.h2h].lact.dropna()
            if len(h) >= 10 and len(t) >= 10:
                R.append(dict(bin=str(b), n_h2h=len(h),
                              lnc_active_h2h=float(h.mean()),
                              n_t2t=len(t), lnc_active_t2t=float(t.mean())))
        B = pd.DataFrame(R)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "I_pair_orientation.tsv")), sep="\t", index=False)
        B.to_csv(safe_write_path(os.path.join(
            OUT, "I_orientation_bins.tsv")), sep="\t", index=False)
        for r in T.itertuples():
            st.note("{}: n={} spacing={:.0f} promoter overlap={:.0f} "
                    "CAGE lnc={:.1%} partner={:.1%}".format(
                        r.subset, r.n, r.median_spacing,
                        r.median_promoter_overlap, r.lnc_cage_active,
                        r.partner_cage_active))
        st.summary_md(os.path.join(
            RES, "step_summaries", "85_pair_orientation.md"),
            "85 divergent pair orientation",
            df_to_md(T, ".4f") + "\n\n" + df_to_md(B, ".4f"))


if __name__ == "__main__":
    main()
