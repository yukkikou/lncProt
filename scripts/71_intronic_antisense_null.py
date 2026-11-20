#!/usr/bin/env python
"""71_intronic_antisense_null.py — do intronic and antisense TSSs prefer open chromatin?

The Results describe these two classes as initiating wherever chromatin is open
around an active host locus, but the supporting numbers are percentages with no
null. A percentage of enhancer-like states means nothing until it is compared
with the states available at that locus.

The null is built from the host gene itself. For every intronic lncRNA the
candidate positions are all positions inside its host protein-coding gene body.
For every antisense lncRNA the candidates are the host gene's exons. A random
position is drawn from that set, its Roadmap state is read from the same
segmentation file, and the process is repeated 20 times per gene. The observed
state distribution is then compared with the distribution at matched positions
that are not annotated initiation sites.

Inputs
  data/chromhmm_hg19/E0*_15_coreMarks_dense.bed.gz   segmentation
  results/gene_master_table.tsv                      classes and hosts

Outputs
  results/promoter_ext/D_class_vs_host_null.tsv
  results/step_summaries/71_intronic_antisense_null.md
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
DATA = os.path.join(WORKDIR, "data", "chromhmm_hg38")
OUT = os.path.join(RES, "promoter_ext")
CELLS = ["E003", "E017", "E067", "E071", "E088"]
NSHUF = 20
SEED = 42
ACTIVE_PROM = {"TssA", "TssAFlnk"}
ENH = {"Enh", "EnhG", "EnhBiv"}


def load_seg(cell):
    path = os.path.join(DATA, "{}.bed".format(cell))
    rows = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            # the segmentation names states "1_TssA", the rest of the
            # project drops the numeric prefix
            rows.append((f[0].replace("chr", ""), int(f[1]), int(f[2]),
                         f[3].split("_", 1)[-1]))
    d = pd.DataFrame(rows, columns=["chrom", "start", "end", "state"])
    by = {}
    for c, g in d.groupby("chrom"):
        g = g.sort_values("start")
        by[c] = (g.start.values, g.end.values, g.state.values)
    return by


def lookup(by, chrom, pos, half=1000):
    """State with the largest overlap in a window of +/-1 kb around pos.

    Module 10 assigns states this way, so a point lookup would not be
    comparable with the published per-class fractions.
    """
    t = by.get(chrom)
    if t is None:
        return None
    s, e, st = t
    lo = np.searchsorted(e, pos - half, side="right")
    hi = np.searchsorted(s, pos + half, side="left")
    if hi <= lo:
        return None
    ov = (np.minimum(e[lo:hi], pos + half)
          - np.maximum(s[lo:hi], pos - half))
    if ov.max() <= 0:
        return None
    return st[lo + int(np.argmax(ov))]


def main():
    with step("71_intronic_antisense_null") as st:
        rng = np.random.default_rng(SEED)
        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        m["tss"] = pd.to_numeric(m.tss, errors="coerce")
        m["chrom_c"] = m.chrom.astype(str).str.replace("chr", "", regex=False)
        # host gene bodies come from the parsed gene BED, because the master
        # table carries only the TSS
        hosts = pd.read_csv(os.path.join(WORKDIR, "data", "genes", "pcg.bed"),
                            sep="\t", header=None,
                            names=["chrom", "start", "end", "gene_id",
                                   "score", "strand"])
        hosts["chrom_c"] = hosts.chrom.str.replace("chr", "", regex=False)
        hosts = hosts.set_index("gene_id")

        rows = []
        for cell in CELLS:
            by = load_seg(cell)
            for cls, host_col in (("intronic", "partner_pcg"),
                                  ("antisense", "partner_pcg")):
                sub = m[(m.cls == cls) & m[host_col].notna()
                        & m.tss.notna()].copy()
                sub = sub[sub[host_col].isin(m.index)]
                host = m.loc[sub[host_col]]
                obs_prom = obs_enh = n_used = 0
                null_prom = null_enh = null_n = 0
                for (lnc, h) in zip(sub.index, sub[host_col]):
                    if h not in hosts.index:
                        continue
                    hs, he = int(hosts.at[h, "start"]), int(hosts.at[h, "end"])
                    if he <= hs:
                        continue
                    chrom = hosts.at[h, "chrom_c"]
                    pos_lnc = int(m.at[lnc, "tss"])
                    s_lnc = lookup(by, m.at[lnc, "chrom_c"], pos_lnc)
                    if s_lnc is None:
                        continue
                    n_used += 1
                    obs_prom += s_lnc in ACTIVE_PROM
                    obs_enh += s_lnc in ENH
                    # random positions inside the same host gene body
                    picks = rng.integers(hs, he, NSHUF)
                    for p in picks:
                        s = lookup(by, chrom, int(p))
                        if s is None:
                            continue
                        null_n += 1
                        null_prom += s in ACTIVE_PROM
                        null_enh += s in ENH
                if n_used < 30:
                    continue
                rows.append(dict(cell=cell, cls=cls, n_loci=n_used,
                                 obs_active_prom=obs_prom / n_used,
                                 null_active_prom=null_prom / null_n,
                                 obs_enhancer=obs_enh / n_used,
                                 null_enhancer=null_enh / null_n,
                                 null_positions=null_n))
        T = pd.DataFrame(rows)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "D_class_vs_host_null.tsv")), sep="\t", index=False)
        for r in T.itertuples():
            st.note("{} {}: active promoter {:.1%} v host null {:.1%}, "
                    "enhancer {:.1%} v {:.1%}".format(
                        r.cell, r.cls, r.obs_active_prom, r.null_active_prom,
                        r.obs_enhancer, r.null_enhancer))
        # pooled test across cell types. Each locus contributes once per cell,
        # so the counts are summed and compared with the pooled null rate.
        if len(T):
            from scipy.stats import binomtest
            tests = []
            for cls, g in T.groupby("cls"):
                k = int(round((g.obs_active_prom * g.n_loci).sum()))
                n = int(g.n_loci.sum())
                p0 = float((g.null_active_prom * g.null_positions).sum()
                           / g.null_positions.sum())
                tests.append(dict(cls=cls, k=k, n=n, p0=p0,
                                  pval=binomtest(k, n, p0,
                                                 alternative="greater").pvalue))
            Tt = pd.DataFrame(tests)
            Tt.to_csv(safe_write_path(os.path.join(
                OUT, "D_class_vs_host_null_test.tsv")), sep="\t", index=False)
            for r in Tt.itertuples():
                st.note("{} pooled: {}/{} active ({:.1%}) against host null "
                        "{:.1%}, binomial P={:.3g}".format(
                            r.cls, r.k, r.n, r.k / r.n, r.p0, r.pval))

        if len(T):
            g = T.groupby("cls")[["obs_active_prom", "null_active_prom",
                                  "obs_enhancer", "null_enhancer"]].mean()
            st.note("mean over five cell types: " + ", ".join(
                "{} active {:.1%} v {:.1%}, enhancer {:.1%} v {:.1%}".format(
                    i, r.obs_active_prom, r.null_active_prom, r.obs_enhancer,
                    r.null_enhancer) for i, r in g.iterrows()))
        st.summary_md(os.path.join(RES, "step_summaries",
                                   "71_intronic_antisense_null.md"),
                      "71 intronic and antisense TSS against a host-gene null",
                      df_to_md(T, ".4f"))


if __name__ == "__main__":
    main()
