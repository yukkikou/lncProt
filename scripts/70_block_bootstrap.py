#!/usr/bin/env python
"""70_block_bootstrap.py — confidence intervals that respect spatial clustering.

LncRNA loci are not spread evenly across the genome. They sit in gene deserts,
in clusters and in segmental duplications, so genes are not independent
observations and an ordinary standard error is too small. Module 45 applied a
chromosome-block bootstrap to the evolutionary statistics only. This module
applies the same resampling to the promoter and activity contrasts, which carry
the largest effects in the paper and therefore the ones least likely to be
artefacts of clustering, and to the tissue specificity contrast, which carries
a small one and is the most likely to be.

Every statistic is a difference of medians between a class and its comparator,
resampled over 10 Mb genomic blocks, 1,000 times, with a fixed seed.

Outputs
  results/evo_boot/bootstrap_promoter.tsv
  results/step_summaries/70_block_bootstrap.md
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "evo_boot")
BLOCK = 10_000_000
NBOOT = 1000
SEED = 42


def add_block(df, chrom="chrom", tss="tss"):
    v = pd.to_numeric(df[tss], errors="coerce")
    df = df.assign(_block=(v // BLOCK).astype("Int64"))
    df["_key"] = df[chrom].astype(str) + ":" + df._block.astype(str)
    return df


def boot_diff(df, a_cls, b_cls, col, rng, agg="median"):
    """Difference of an aggregate, resampled over genomic blocks.

    A proportion needs the mean; a continuous trait needs the median.
    """
    d = df[[col, "cls", "_key"]].dropna()
    a = d[d.cls == a_cls]
    b = d[d.cls == b_cls]
    if len(a) < 50 or len(b) < 50:
        return None
    f = np.mean if agg == "mean" else np.median
    obs = float(f(a[col]) - f(b[col]))
    keys = np.array(sorted(set(d._key)))
    ia = a.groupby("_key")[col].apply(list).to_dict()
    ib = b.groupby("_key")[col].apply(list).to_dict()
    vals = []
    for _ in range(NBOOT):
        k = rng.choice(keys, len(keys), replace=True)
        va = [x for kk in k for x in ia.get(kk, ())]
        vb = [x for kk in k for x in ib.get(kk, ())]
        if len(va) > 10 and len(vb) > 10:
            vals.append(f(va) - f(vb))
    if len(vals) < 100:
        return None
    return dict(comparison="{} minus {}".format(a_cls, b_cls), trait=col,
                agg=agg,
                n_a=len(a), n_b=len(b), observed=obs,
                ci_low=float(np.percentile(vals, 2.5)),
                ci_high=float(np.percentile(vals, 97.5)),
                n_boot=len(vals), blocks=len(keys))


def main():
    with step("70_block_bootstrap") as st:
        rng = np.random.default_rng(SEED)
        seq = pd.read_csv(os.path.join(RES, "seq_features.tsv"), sep="\t",
                          low_memory=False)
        seq = add_block(seq)
        rows = []
        for col, agg in (("gc", "median"), ("cpg_oe", "median"),
                         ("cpg_island", "mean")):
            if col not in seq:
                continue
            for a, b in (("divergent", "intergenic"), ("intergenic", "PCG"),
                         ("divergent", "PCG")):
                r = boot_diff(seq, a, b, col, rng, agg)
                if r:
                    rows.append(r)
        # activity
        c = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                        low_memory=False)
        c = add_block(c)
        c["active"] = c.cage_refined.fillna(0).astype(float)
        for a, b in (("divergent", "intergenic"), ("divergent", "PCG")):
            r = boot_diff(c, a, b, "active", rng, "mean")
            if r:
                rows.append(r)
        # tissue specificity, the small effect
        g = pd.read_csv(os.path.join(
            RES, "tissue_disease", "tables", "gtex_gene_tau.tsv"), sep="\t")
        g = g.rename(columns={"Unnamed: 0": "gene_id"}).dropna(
            subset=["tau_gtex"])
        # the GTEx table carries version-free identifiers, the sequence table
        # versioned ones, so the join needs both stripped
        seq_base = seq.assign(base=seq.gene_id.str.split(".").str[0]
                              ).drop_duplicates("base")[["base", "chrom", "tss"]]
        g["base"] = g.gene_id.str.split(".").str[0]
        g = g.merge(seq_base, on="base", how="left")
        g = add_block(g)
        for a, b in (("divergent", "PCG"), ("intergenic", "PCG")):
            r = boot_diff(g, a, b, "tau_gtex", rng)
            if r:
                rows.append(r)

        T = pd.DataFrame(rows)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "bootstrap_promoter.tsv")), sep="\t", index=False)
        for r in T.itertuples():
            st.note("{} {}: {:.3f} [{:.3f}, {:.3f}] ({} blocks)".format(
                r.trait, r.comparison, r.observed, r.ci_low, r.ci_high,
                r.blocks))
        st.summary_md(os.path.join(RES, "step_summaries",
                                   "70_block_bootstrap.md"),
                      "70 block bootstrap for promoter and activity contrasts",
                      df_to_md(T, ".3f"))


if __name__ == "__main__":
    main()
