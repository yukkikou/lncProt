#!/usr/bin/env python
"""12a_pub_stats.py — publication-grade statistical tests, class vs PCG.

Continuous features: Mann-Whitney U + rank-biserial effect size (Cliff's delta)
Binary features:     Fisher exact + log2OR
All p-values BH-FDR within each feature family.

Inputs : results/gene_master_table.tsv (joined per-gene features are read from
         the individual module tables to avoid column drift)
Outputs: results/pub_stats.tsv, step_summaries/12a_pub_stats.md
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic", "pseudogene", "transcribed_pseudogene"]


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def build_frame():
    m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                    index_col=0, low_memory=False)
    for f, cols in [("seq_features.tsv", ["gc", "cpg_oe", "TATA", "Inr",
                                          "CCAAT", "prom_class", "cpg_island"]),
                    ("te_promoter_stats.tsv", ["te_cov", "tss_in_te", "tss_in_ltr"]),
                    ("conservation.tsv", ["phast_mean"]),
                    ("cage_activity.tsv", ["active", "max_tpm"])]:
        d = pd.read_csv(os.path.join(RES, f), sep="\t", index_col=0,
                        low_memory=False)
        keep = [c for c in cols if c in d.columns and c not in m.columns]
        if keep:
            m = m.join(d[keep], how="left")
    # chromhmm prom-like (mean across cells) as binary-ish
    ch = pd.read_csv(os.path.join(RES, "chromhmm_states.tsv"), sep="\t",
                     index_col=0, low_memory=False)
    m["prom_like_frac"] = ch.eq("TssA").mean(axis=1).reindex(m.index).values
    m["log10_tpm"] = np.log10(m["max_tpm"].astype(float).fillna(0) + 1)
    return m


def main():
    with step("12a_pub_stats") as st:
        m = build_frame()
        pcg = m[m.group == "pcg"]

        CONT = {"gc": "GC content", "cpg_oe": "CpG obs/exp",
                "phast_mean": "phastCons (470w)", "te_cov": "TE coverage",
                "log10_tpm": "log10 max TPM (promoter output)"}
        BIN = {"TATA": "TATA box", "Inr": "Initiator", "CCAAT": "CCAAT box",
               "cpg_island": "CpG island overlap", "HCP": "HCP promoter",
               "tss_in_ltr": "TSS inside LTR/ERV", "tss_in_te": "TSS inside any TE",
               "active": "CAGE-active promoter", "enhancer_overlap": "FANTOM enhancer at promoter"}

        rows = []
        # continuous
        for feat, label in CONT.items():
            x = pcg[feat].dropna().values
            for cls in CLASSES:
                sub = m[m.cls == cls]
                y = sub[feat].dropna().values
                if len(y) < 20:
                    continue
                u, p = mannwhitneyu(y, x, alternative="two-sided")
                delta = 2 * u / (len(y) * len(x)) - 1     # Cliff's delta
                rows.append(dict(feature=label, type="continuous", cls=cls,
                                 n_class=len(y),
                                 mean_class=float(np.mean(y)), mean_PCG=float(np.mean(x)),
                                 cliffs_delta=float(delta), pval=p))
        # binary
        for feat, label in BIN.items():
            if feat == "HCP":
                pc = (pcg["prom_class"] == "HCP")
                cl_vals = (m["prom_class"] == "HCP")
            else:
                pc = pcg[feat].fillna(0).astype(bool)
                cl_vals = m[feat].fillna(0).astype(bool)
            n11_ref = int(pc.sum())
            n01_ref = int((~pc).sum())
            for cls in CLASSES:
                sub_mask = (m.cls == cls)
                a = (cl_vals & sub_mask).sum()
                b = (sub_mask & ~cl_vals).sum()
                if a + b < 20:
                    continue
                orr, p = fisher_exact([[int(a), int(b)], [n11_ref, n01_ref]])
                rows.append(dict(feature=label, type="binary", cls=cls,
                                 n_class=int(a + b),
                                 mean_class=float(a / max(a + b, 1)),
                                 mean_PCG=float(n11_ref / max(n11_ref + n01_ref, 1)),
                                 log2OR=float(np.log2(orr)) if orr > 0 else np.nan,
                                 pval=p))
        t = pd.DataFrame(rows)
        t["fdr"] = bh_fdr(t["pval"].values)
        t = t.sort_values(["feature", "cls"])
        t.to_csv(safe_write_path(os.path.join(RES, "pub_stats.tsv")),
                 sep="\t", index=False)
        st.note(f"tests={len(t)}; significant(FDR<0.05)={int((t.fdr < 0.05).sum())}")

        # compact view for the paper: key features x classes
        key_feats = ["GC content", "CpG obs/exp", "phastCons (470w)",
                     "TSS inside LTR/ERV", "TATA box", "CAGE-active promoter",
                     "TE coverage"]
        pv = t[t.feature.isin(key_feats)]
        sig = pv.assign(star=lambda d: np.where(d.fdr < 1e-10, "***",
                        np.where(d.fdr < 1e-3, "**",
                        np.where(d.fdr < 0.05, "*", "ns"))))
        body = ("## Class vs PCG statistical tests (all classes, key features)\n\n"
                + df_to_md(sig, ".3f", max_rows=60)
                + "\n\n Cliff's delta: |d|>0.33 medium, >0.47 large. "
                  "log2OR for binary features. All p BH-FDR corrected.")
        st.summary_md("results/step_summaries/12a_pub_stats.md",
                      "12a publication statistics", body)
        # quick highlights
        for cls in ("divergent", "intergenic"):
            r = t[(t.cls == cls) & (t.fdr < 0.05)]
            st.note(f"{cls}: {len(r)} significant tests of {len(t[t.cls==cls])}")


if __name__ == "__main__":
    main()
