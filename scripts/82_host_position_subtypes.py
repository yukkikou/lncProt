#!/usr/bin/env python
"""82_host_position_subtypes.py — where inside the host locus do the middle
classes initiate?

Intronic and antisense lncRNAs may be promoter-proximal read-through rather
than opportunistic initiation. This module stratifies them by position
relative to the host coding gene and compares CAGE support, TSS-in-LTR and
conservation across strata.

Intronic strata are promoter-proximal (within 2 kb of host TSS), first 10% of
the gene body, and distal. Antisense strata are 5-prime (within 2 kb of host
TSS), internal and 3-prime (within 2 kb of host end).

Outputs
  results/revision_checks/F_host_position.tsv
  results/revision_checks/F_host_position_tests.tsv
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


def main():
    with step("82_host_position_subtypes") as st:
        os.makedirs(OUT, exist_ok=True)
        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        m["tss"] = pd.to_numeric(m.tss, errors="coerce")
        hosts = pd.read_csv(os.path.join(WORKDIR, "data", "genes", "pcg.bed"),
                            sep="\t", header=None,
                            names=["chrom", "start", "end", "gene_id",
                                   "score", "strand"]).set_index("gene_id")
        cage = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                           low_memory=False)
        cage["base"] = cage.gene_id.astype(str).str.split(".").str[0]
        cage = cage.drop_duplicates("base").set_index("base")
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                         index_col=0, low_memory=False)[["tss_in_ltr"]]
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)[["phast_mean"]]
        m["base"] = m.index.astype(str).str.split(".").str[0]
        m["active"] = m.base.map(cage.active)
        m = m.join(te).join(ph)

        rows = []
        for cls in ("intronic", "antisense"):
            sub = m[(m.cls == cls) & m.partner_pcg.notna()
                    & m.tss.notna()].copy()
            sub["host_base"] = sub.partner_pcg.astype(str)
            sub = sub[sub.host_base.isin(hosts.index)]
            hs = hosts.reindex(sub.host_base)
            sub["hstart"] = hs.start.values.astype(float)
            sub["hend"] = hs.end.values.astype(float)
            sub["hstrand"] = hs.strand.values
            sub["htss"] = np.where(sub.hstrand == "+", sub.hstart,
                                   sub.hend)
            body = (sub.hend - sub.hstart).clip(lower=1)
            dist = np.where(sub.hstrand == "+", sub.tss - sub.hstart,
                            sub.hend - sub.tss)
            sub["frac_body"] = np.clip(dist / body, 0, 1)
            if cls == "intronic":
                sub["stratum"] = np.where(
                    (sub.tss - sub.htss).abs() <= 2000, "promoter_proximal",
                    np.where(sub.frac_body <= 0.10, "first_10pct",
                             "distal"))
            else:
                sub["stratum"] = np.select(
                    [(sub.tss - sub.htss).abs() <= 2000,
                     ((sub.tss - sub.hend).abs() <= 2000)
                     | ((sub.tss - sub.hstart).abs() <= 2000)],
                    ["5prime", "3prime"], default="internal")
            for stratum, g in sub.groupby("stratum"):
                rows.append(dict(cls=cls, stratum=stratum, n=len(g),
                                 frac_cage_active=float(g.active.mean()),
                                 frac_tss_in_ltr=float(g.tss_in_ltr.mean()),
                                 median_phast=float(g.phast_mean.median()),
                                 median_frac_body=float(
                                     g.frac_body.median())))
        T = pd.DataFrame(rows)

        tests = []
        for cls in ("intronic", "antisense"):
            sub = T[T.cls == cls]
            if len(sub) < 2:
                continue
            a = sub.iloc[0]
            b = sub.iloc[-1]
            if a.n > 20 and b.n > 20:
                tests.append(dict(cls=cls, contrast="{} vs {}".format(
                    a.stratum, b.stratum),
                    frac_cage_active_diff=a.frac_cage_active
                    - b.frac_cage_active,
                    frac_tss_in_ltr_diff=a.frac_tss_in_ltr
                    - b.frac_tss_in_ltr,
                    phast_diff=a.median_phast - b.median_phast))
        U = pd.DataFrame(tests)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "F_host_position.tsv")), sep="\t", index=False)
        if len(U):
            U.to_csv(safe_write_path(os.path.join(
                OUT, "F_host_position_tests.tsv")), sep="\t", index=False)
        for r in T.itertuples():
            st.note("{} {}: n={} CAGE={:.1%} LTR={:.1%} phast={:.3f}".format(
                r.cls, r.stratum, r.n, r.frac_cage_active,
                r.frac_tss_in_ltr, r.median_phast))
        st.summary_md(os.path.join(
            RES, "step_summaries", "82_host_position_subtypes.md"),
            "82 host-position subtypes",
            df_to_md(T, ".4f") + ("\n\n" + df_to_md(U, ".4f")
                                  if len(U) else ""))


if __name__ == "__main__":
    main()
