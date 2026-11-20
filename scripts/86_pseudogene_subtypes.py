#!/usr/bin/env python
"""86_pseudogene_subtypes.py — split the pseudogene comparator by origin.

The manuscript pools 15,210 pseudogenes. Processed pseudogenes arise by
retrotransposition and have no ancestral promoter. Duplicated pseudogenes
retain one. Pooling them guarantees a compositionally extreme mixture, so the
decayed-syntax conclusion is partly built into the group definition. This
module separates the two origins, and transcribed from non-transcribed loci,
and re-runs the promoter comparison for each subset.

Outputs
  results/revision_checks/J_pseudogene_subtypes.tsv
  results/revision_checks/J_pseudogene_tests.tsv
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "revision_checks")


def origin(gt):
    gt = str(gt).lower()
    if "processed" in gt:
        return "processed"
    if "unprocessed" in gt or "unitary" in gt:
        return "duplicated"
    return "other"


def main():
    with step("86_pseudogene_subtypes") as st:
        os.makedirs(OUT, exist_ok=True)
        genes = pd.read_csv(os.path.join(WORKDIR, "data", "genes", "genes.tsv"),
                            sep="\t").set_index("gene_id")
        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        m["base"] = m.index.astype(str).str.split(".").str[0]
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                         index_col=0, low_memory=False)
        sf = pd.read_csv(os.path.join(RES, "seq_features.tsv"), sep="\t",
                         low_memory=False)
        sf["base"] = sf.gene_id.astype(str).str.split(".").str[0]
        sf = sf.drop_duplicates("base").set_index("base")[["gc", "cpg_oe"]]
        cage = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                           low_memory=False)
        cage["base"] = cage.gene_id.astype(str).str.split(".").str[0]
        cage = cage.drop_duplicates("base").set_index("base")[["active"]]
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)
        ph["base"] = ph.index.astype(str).str.split(".").str[0]
        ph = ph.drop_duplicates("base").set_index("base")[["phast_mean"]]

        rows = []
        for cls in ("pseudogene", "transcribed_pseudogene"):
            sub = m[m.cls == cls].copy()
            sub = sub[sub.partner_pcg.notna() | True]
            sub["origin"] = [origin(genes.gene_type.get(g, ""))
                             for g in sub.index]
            for (org, trans), g in sub.groupby(
                    ["origin", "transcribed_pseudo"]):
                if len(g) < 20:
                    continue
                lab = "{} {}".format(
                    "transcribed" if trans == 1 else "untranscribed", org)
                rows.append(dict(
                    cls=cls, subset=lab, n=len(g),
                    frac_cage_active=float(
                        g.base.map(cage.active).mean()),
                    frac_high_cpg=float(
                        (g.base.map(sf.cpg_oe) >= 0.6).mean()),
                    mean_gc=float(g.base.map(sf.gc).mean()),
                    mean_te_cov=float(
                        te.reindex(g.index).te_cov.mean()),
                    frac_tss_in_te=float(
                        te.reindex(g.index).tss_in_te.mean()),
                    frac_tss_in_ltr=float(
                        te.reindex(g.index).tss_in_ltr.mean()),
                    median_phast=float(
                        g.base.map(ph.phast_mean).median())))
        T = pd.DataFrame(rows)

        # tests against PCG background
        pcg = m[m.group == "pcg"]
        bg = dict(
            cage=float(pcg.base.map(cage.active).mean()),
            hcp=float((pcg.base.map(sf.cpg_oe) >= 0.6).mean()),
            te=float(te.reindex(pcg.index).te_cov.mean()),
            tltr=float(te.reindex(pcg.index).tss_in_ltr.mean()),
            phast=float(pcg.base.map(ph.phast_mean).median()))
        tests = []
        for r in T.itertuples():
            sub = m[m.cls == r.cls][m[m.cls == r.cls].base.isin(
                m[m.cls == r.cls].base)]
            tests.append(dict(subset=r.subset, n=r.n,
                              or_cage_vs_pcg=float(
                                  fisher_exact([[round(r.frac_cage_active * r.n),
                                                 r.n - round(r.frac_cage_active * r.n)],
                                                [round(bg["cage"] * len(pcg)),
                                                 len(pcg) - round(bg["cage"] * len(pcg))]])[0]),
                              te_cov_vs_pcg=r.mean_te_cov - bg["te"],
                              phast_vs_pcg=r.median_phast - bg["phast"]))
        S = pd.DataFrame(tests)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "J_pseudogene_subtypes.tsv")), sep="\t", index=False)
        S.to_csv(safe_write_path(os.path.join(
            OUT, "J_pseudogene_tests.tsv")), sep="\t", index=False)
        st.note("PCG background: CAGE {:.1%}, TE {:.3f}, phast {:.3f}".format(
            bg["cage"], bg["te"], bg["phast"]))
        for r in T.itertuples():
            st.note("{}: n={} CAGE={:.1%} TE={:.3f} LTR={:.1%} 5mC={:.3f}"
                    .format(r.subset, r.n, r.frac_cage_active, r.mean_te_cov,
                            r.frac_tss_in_ltr, r.median_phast))
        st.summary_md(os.path.join(
            RES, "step_summaries", "86_pseudogene_subtypes.md"),
            "86 pseudogene subtypes",
            df_to_md(T, ".4f") + "\n\n" + df_to_md(S, ".4f"))


if __name__ == "__main__":
    main()
