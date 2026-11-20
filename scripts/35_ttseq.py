#!/usr/bin/env python
"""35_ttseq.py — initiation vs stability by birth class using TT-seq.

GSE75792 (Schwalb et al.) transcriptional units (TUs) carry per-TU synthesis
and decay rates. We map each TU TSS to the nearest GENCODE v46 gene TSS
(same strand, <=500 bp), assign birth class, and compare synthesis (initiation
proxy), decay (stability proxy) and their ratio between classes.
Outputs -> results/ttseq/*.tsv + step summary.
"""
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "ttseq")
GTF = os.path.join(WORKDIR, "data", "perturb", "ttseq",
                   "GSE75792_transcript.annotation.gtf.gz")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def mwu(y, x):
    u, p = mannwhitneyu(y, x, alternative="two-sided")
    return float(2 * u / (len(y) * len(x)) - 1), float(p)


def parse_tu():
    rows = []
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            attr = p[8]
            ms = re.search(r'synthesis rate "([0-9.eE+-]+)"', attr)
            md = re.search(r'decay rate "([0-9.eE+-]+)"', attr)
            if not ms or not md:
                continue
            rows.append(dict(chrom=p[0], start=int(p[3]), end=int(p[4]),
                             strand=p[6], feature=p[2],
                             syn=float(ms.group(1)),
                             decay=float(md.group(1)),
                             gencode='GENCODE "TRUE"' in attr))
    return pd.DataFrame(rows)


def main():
    with step("35_ttseq") as st:
        tu = parse_tu()
        tu["tss"] = np.where(tu.strand == "+", tu.start, tu.end)
        tu["len"] = tu.end - tu.start
        st.note("TUs parsed: {} (GENCODE TRUE={})".format(
            len(tu), int(tu.gencode.sum())))

        g = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                        sep="\t", index_col=0, low_memory=False)
        g["tss"] = g["tss"].astype(int)
        g = g[["chrom", "tss", "strand", "cls", "group"]]
        by = {}
        for (c, s), sub in g.groupby(["chrom", "strand"]):
            by[(c, s)] = np.sort(sub.tss.values)
        # map TU tss -> nearest gene tss (same chrom/strand, <=500 bp)
        tu["tss_gene"] = np.nan
        for (c, s), sub in tu.groupby(["chrom", "strand"]):
            arr = by.get((c, s))
            if arr is None:
                continue
            pos = sub.tss.values
            i = np.searchsorted(arr, pos)
            best = np.full(len(pos), np.nan)
            for k in range(len(pos)):
                cand = []
                if i[k] < len(arr):
                    cand.append(arr[i[k]])
                if i[k] > 0:
                    cand.append(arr[i[k] - 1])
                if cand:
                    j = int(np.argmin([abs(x - pos[k]) for x in cand]))
                    d = abs(cand[j] - pos[k])
                    if d <= 500:
                        best[k] = cand[j]
            tu.loc[sub.index, "tss_gene"] = best
        # map tss_gene -> gene id (unique per chrom/strand/tss; keep first)
        key = g.copy()
        key["k"] = key.chrom + ":" + key.strand + ":" + key.tss.astype(str)
        kmap = key.drop_duplicates("k").set_index("k")
        tu["key"] = (tu.chrom + ":" + tu.strand + ":" +
                     pd.Series(tu.tss_gene).astype("Int64").astype(str))
        tu = tu.join(kmap[["cls", "group"]], on="key")
        mapped = tu.dropna(subset=["cls"]).copy()
        # one TU per gene: closest TSS; keep longest on tie
        mapped["dist"] = (mapped.tss - mapped.tss_gene).abs()
        mapped = mapped.sort_values(["key", "dist", "len"],
                                    ascending=[True, True, False])
        per = mapped.drop_duplicates("key").copy()
        g2 = g.reset_index()
        g2["k"] = g2.chrom + ":" + g2.strand + ":" + g2.tss.astype(str)
        gid = g2.drop_duplicates("k").set_index("k")["gene_id"]
        per["gene_id"] = per["key"].map(gid)
        per["log_syn"] = np.log10(per.syn + 1e-6)
        per["log_decay"] = np.log10(per.decay + 1e-6)
        per["log_ratio"] = per.log_syn - per.log_decay
        per["syn_per_kb"] = per.syn / (per["len"].clip(lower=200) / 1000.0)
        per.to_csv(safe_write_path(os.path.join(
            OUT, "TU_mapped.tsv")), sep="\t", index=False)
        st.note("mapped TUs: {} ({:.1%}); genes: {}".format(
            len(per), len(per) / len(tu), per.gene_id.nunique()))

        rows = []
        pcg = per[per.group == "pcg"]
        for feat in ["log_syn", "log_decay", "log_ratio", "syn_per_kb"]:
            x = pcg[feat].dropna().values
            for cls in CLASSES + ["pseudogene", "transcribed_pseudogene"]:
                y = per.loc[per.cls == cls, feat].dropna().values
                if len(y) < 20:
                    continue
                d, p = mwu(y, x)
                rows.append(dict(feature=feat, cls=cls, n=len(y),
                                 median_class=float(np.median(y)),
                                 median_PCG=float(np.median(x)),
                                 cliffs_delta=d, pval=p))
        T = pd.DataFrame(rows)
        T["fdr"] = bh_fdr(T.pval.values)
        T.to_csv(safe_write_path(os.path.join(OUT, "tests.tsv")),
                 sep="\t", index=False)

        # intergenic LTR vs non-LTR
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                         index_col=0, low_memory=False)[["tss_in_ltr"]]
        per2 = per.set_index("gene_id").join(te)
        ig = per2[per2.cls == "intergenic"]
        lrows = []
        for feat in ["log_syn", "log_decay", "log_ratio"]:
            a = ig.loc[ig.tss_in_ltr == 1, feat].dropna().values
            b = ig.loc[ig.tss_in_ltr != 1, feat].dropna().values
            if len(a) >= 20 and len(b) >= 20:
                d, p = mwu(a, b)
                lrows.append(dict(feature=feat, n=len(a) + len(b),
                                  median_LTR=float(np.median(a)),
                                  median_nonLTR=float(np.median(b)),
                                  cliffs_delta=d, pval=p))
        L = pd.DataFrame(lrows)
        L.to_csv(safe_write_path(os.path.join(OUT, "intergenic_LTR.tsv")),
                 sep="\t", index=False)

        act = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                          index_col=0, low_memory=False)[
            ["max_tpm", "tau", "active"]]
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)[
            ["phylop_mean"]]
        corr = per.set_index("gene_id").join(act).join(ph)
        crows = []
        for feat in ["log_syn", "log_ratio"]:
            for other in ["max_tpm", "tau", "phylop_mean"]:
                ok = corr[[feat, other]].dropna()
                if len(ok) > 50:
                    r, p = spearmanr(ok[feat], ok[other])
                    crows.append(dict(feature=feat, other=other, n=len(ok),
                                      rho=float(r), pval=float(p)))
        C = pd.DataFrame(crows)
        st.note("TT-seq: mapped genes={}; " + "; ".join(
            "{} {} d={:.3f} p={:.1e}".format(
                r.feature, r.cls, r.cliffs_delta, r.pval)
            for r in T.itertuples() if r.cls in
            ("divergent", "intergenic")))
        st.summary_md("results/ttseq/step_summaries/35_ttseq.md",
                      "35 TT-seq initiation vs stability",
                      "## Class vs PCG (per mapped TU)\n\n"
                      + df_to_md(T, ".3f", max_rows=40)
                      + "\n\n## Intergenic LTR vs non-LTR\n\n"
                      + df_to_md(L, ".3f")
                      + "\n\n## Correlations\n\n" + df_to_md(C, ".3f")
                      + "\n\nSynthesis = initiation proxy; decay = stability "
                        "proxy; log_ratio = synthesis - decay.")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    main()
