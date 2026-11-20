#!/usr/bin/env python
"""79_anage_cancer_rates.py — the lifespan axis at lncRNA-associated loci.

Per-species alignments are not on this cluster, so the 58-species rate
analysis of Zhang et al. cannot be recomputed from raw data. The paper's
Source Data do provide per-site evolutionary rates at cancer-associated
positions, separately for their longer-lived and shorter-lived species groups,
keyed by Ensembl transcript. This module aggregates those rates per coding
gene and asks whether the long-lived conservation signal is stronger at the
genes that lncRNAs borrow (divergent partners) or sit inside (intronic and
antisense hosts), and whether it concentrates in Cancer Gene Census genes.

Outputs
  results/comparative/anage_gene_rates.tsv
  results/comparative/anage_rate_tests.tsv
  results/step_summaries/79_anage_cancer_rates.md
"""
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
SRC = os.path.join(WORKDIR, "data", "anage", "paper_source", "extracted")
OUT = os.path.join(RES, "comparative")
GTF = ("/media/share/resource/reference/homo_sapiens/GRCh38/annotation/"
       "gencode.v46.annotation.gtf")
LONG = os.path.join(SRC, "Fig1G-R and Fig 2 BoxPlotData_long.cancer.rate.txt")
SHORT = os.path.join(SRC, "Fig1G-R and Fig 2 BoxPlotData_short.cancer.rate.txt")
MIN_SITES = 5


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def transcript_map():
    rows = []
    for line in open(GTF):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "transcript":
            continue
        tid = re.search(r'transcript_id "([^"]+)"', f[8])
        gid = re.search(r'gene_id "([^"]+)"', f[8])
        if tid and gid:
            rows.append((tid.group(1).split(".")[0],
                         gid.group(1).split(".")[0]))
    return pd.DataFrame(rows, columns=["enst", "ensg"]).drop_duplicates("enst")


def main():
    with step("79_anage_cancer_rates") as st:
        long = pd.read_csv(LONG, sep="\t", header=None,
                           names=["enst", "pos", "rate"], dtype={0: str})
        short = pd.read_csv(SHORT, sep="\t", header=None,
                            names=["enst", "pos", "rate"], dtype={0: str})
        st.note("rate rows long={:,} short={:,}".format(len(long), len(short)))
        d = long.merge(short, on=["enst", "pos"], suffixes=("_long", "_short"))
        st.note("matched sites={:,}".format(len(d)))

        tm = transcript_map()
        d = d.merge(tm, on="enst", how="left").dropna(subset=["ensg"])
        G = d.groupby("ensg").agg(
            n_sites=("rate_long", "size"),
            rate_long=("rate_long", "mean"),
            rate_short=("rate_short", "mean"))
        G = G[G.n_sites >= MIN_SITES].copy()
        G["delta_short_minus_long"] = G.rate_short - G.rate_long
        st.note("genes with >= {} sites={}".format(MIN_SITES, len(G)))

        # annotation
        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        m["base"] = m.index.astype(str).str.split(".").str[0]
        pcg = m[m.group == "pcg"].copy().set_index("base")
        cg = pd.read_csv(os.path.join(
            WORKDIR, "data", "cosmic", "extracted",
            "Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz"), sep="\t",
            low_memory=False)
        cg["ens"] = cg.SYNONYMS.astype(str).str.extract(r"(ENSG\d{11})")[0]
        role = cg.ROLE_IN_CANCER.fillna("").astype(str)
        for lab, mask in (("cgc", pd.Series(True, index=cg.index)),
                          ("onc", role.str.contains("oncogene", case=False)),
                          ("tsg", role.str.contains("TSG", case=False)),
                          ("fusion", role.str.contains("fusion", case=False))):
            sy = set(cg.loc[mask, "GENE_SYMBOL"].astype(str))
            en = set(cg.loc[mask, "ens"].dropna())
            pcg[lab] = pcg.gene_name.astype(str).isin(sy) | pcg.index.isin(en)

        dp = pd.read_csv(os.path.join(RES, "divergent_pairs.tsv"), sep="\t")
        partnered = set(dp.pcg_gid.astype(str).str.split(".").str[0])
        hosts = set()
        for cls in ("intronic", "antisense"):
            hosts |= set(m[m.cls == cls].partner_pcg.dropna()
                         .astype(str).str.split(".").str[0])
        pcg["partner"] = pcg.index.isin(partnered)
        pcg["host"] = pcg.index.isin(hosts)

        J = G.join(pcg[["cgc", "onc", "tsg", "fusion", "partner", "host"]],
                   how="inner")
        J = J.rename(columns={"cgc": "is_cgc", "onc": "is_oncogene",
                              "tsg": "is_tsg", "fusion": "is_fusion_role",
                              "partner": "is_divergent_partner",
                              "host": "is_intronic_antisense_host"})
        J.to_csv(safe_write_path(os.path.join(
            OUT, "anage_gene_rates.tsv")), sep="\t", index_label="gene_id")

        rows = []

        def add(test, subset, stat, p, n, note):
            rows.append(dict(test=test, subset=subset, statistic=stat,
                             pval=p, n=n, note=note))

        def cliff(a, b):
            a = np.asarray(a); b = np.asarray(b)
            if len(a) < 5 or len(b) < 5:
                return np.nan
            from scipy.stats import mannwhitneyu as mwu
            u, _ = mwu(a, b)
            return 2 * u / (len(a) * len(b)) - 1

        # overall replication
        w, p = wilcoxon(J.delta_short_minus_long)
        add("lifespan", "short minus long rate, all genes",
            float(J.delta_short_minus_long.median()), p, len(J),
            "positive means higher rate (less conservation) in short-lived")

        for lab, labname in (
                ("is_cgc", "CGC"), ("is_oncogene", "oncogene"),
                ("is_tsg", "TSG"), ("is_fusion_role", "fusion role"),
                ("is_divergent_partner", "divergent partner"),
                ("is_intronic_antisense_host", "intronic/antisense host")):
            a = J[J[lab]].delta_short_minus_long
            b = J[~J[lab]].delta_short_minus_long
            if len(a) > 20 and len(b) > 20:
                u, p = mannwhitneyu(a, b)
                add("lifespan", "delta by {}".format(labname),
                    cliff(a, b), p, len(a) + len(b),
                    "median {:.4f} v {:.4f}; Cliff's delta".format(
                        a.median(), b.median()))

        T = pd.DataFrame(rows)
        T["fdr"] = bh_fdr(T.pval.values)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "anage_rate_tests.tsv")), sep="\t", index=False)
        for r in T.itertuples():
            st.note("{} | {} | stat={:.3f} n={} FDR={:.3g} | {}".format(
                r.test, r.subset, r.statistic, r.n, r.fdr, r.note))
        st.summary_md(os.path.join(
            RES, "step_summaries", "79_anage_cancer_rates.md"),
            "79 AnAge cancer-site rates at lncRNA-associated genes",
            df_to_md(T, ".4f"))


if __name__ == "__main__":
    main()
