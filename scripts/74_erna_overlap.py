#!/usr/bin/env python
"""74_erna_overlap.py — do class TSSs sit in transcribed enhancers (eRNAs)?

The manuscript argues that intronic and antisense lncRNAs initiate where the
host locus is already open, and that intergenic lncRNAs often use retroviral
promoters. eRNAbase, available locally, maps enhancer RNAs across human cell
types and tissues. If the open-host-locus model is right, those two classes
should overlap eRNAs more often than their own neighbourhood predicts, not just
more often than coding promoters.

Two nulls, because eRNA density is spatially clustered:

  local  a random position within +/- 100 kb of the same TSS. This preserves
         the local eRNA environment, so a class that simply lives in eRNA-rich
         regions cannot pass.
  host   a random position inside the same host PCG body, for intronic and
         antisense lncRNAs only. This is the null used in module 71.

Each locus is compared with its own 20 shuffled positions, and the pooled test
is the exact Poisson-binomial over genes (as in module 65), so neither eRNA
density nor gene length drives the result.

Inputs
  /media/share/resource/database/eRNA_annotation/eRNAbase_eRNA_merge_human_nonred.bed
  results/gene_master_table.tsv
  data/genes/pcg.bed
  results/te_promoter_stats.tsv

Outputs
  results/promoter_ext/E_erna_overlap.tsv
  results/step_summaries/74_erna_overlap.md
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.signal import fftconvolve
from scipy.stats import binom

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
ERNA = ("/media/share/resource/database/eRNA_annotation/"
        "eRNAbase_eRNA_merge_human_nonred.bed")
OUT = os.path.join(RES, "promoter_ext")
ORDER = ["PCG", "divergent", "antisense", "intronic", "sense_overlapping",
         "intergenic", "pseudogene", "transcribed_pseudogene"]
NSHUF = 20
HALF = 500
SEED = 42


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def load_intervals():
    """Merged, sorted eRNA intervals per chromosome."""
    d = pd.read_csv(ERNA, sep="\t", header=None, usecols=[0, 1, 2],
                    names=["chrom", "start", "end"])
    d["chrom"] = d.chrom.astype(str).str.replace("chr", "", regex=False)
    out = {}
    for c, g in d.groupby("chrom"):
        s = g.start.to_numpy(dtype=np.int64)
        e = g.end.to_numpy(dtype=np.int64)
        o = np.argsort(s)
        out[c] = (s[o], e[o])
    return out


def overlap(intervals, chrom, pos, half=HALF):
    """True if the window around pos touches any merged eRNA interval."""
    t = intervals.get(chrom)
    if t is None:
        return False
    s, e = t
    j = int(np.searchsorted(s, pos + half, side="right")) - 1
    if j < 0:
        return False
    return bool(e[j] >= pos - half)


def poisson_binomial_p(k_counts, obs, n_shuf):
    """Exact P(T >= obs) when each gene is exchangeable with its shuffles."""
    pmf = np.array([1.0])
    for k in range(1, n_shuf + 2):
        c = int(k_counts[k])
        if c == 0:
            continue
        pmf = fftconvolve(pmf, binom.pmf(np.arange(c + 1), c,
                                         k / (n_shuf + 1.0)))
    pmf = np.clip(pmf, 0, None)
    s = pmf.sum()
    if s > 0:
        pmf /= s
    if obs >= len(pmf):
        return 0.0
    return float(min(max(pmf[obs:].sum(), 0.0), 1.0))


def main():
    with step("74_erna_overlap") as st:
        rng = np.random.default_rng(SEED)
        iv = load_intervals()
        st.note("eRNA intervals: {} chromosomes".format(len(iv)))

        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        m["tss"] = pd.to_numeric(m.tss, errors="coerce")
        m["chrom_c"] = m.chrom.astype(str).str.replace("chr", "", regex=False)
        m = m.dropna(subset=["tss"])
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                         index_col=0, low_memory=False)[["tss_in_ltr"]]
        m = m.join(te)
        hosts = pd.read_csv(os.path.join(WORKDIR, "data", "genes", "pcg.bed"),
                            sep="\t", header=None,
                            names=["chrom", "start", "end", "gene_id",
                                   "score", "strand"])
        hosts["chrom_c"] = hosts.chrom.str.replace("chr", "", regex=False)
        hosts = hosts.set_index("gene_id")

        rows = []
        for cls in ORDER:
            sub = m[m.cls == cls]
            if len(sub) < 20:
                continue
            obs = np.array([overlap(iv, r.chrom_c, int(r.tss))
                            for r in sub.itertuples()])
            # local null: random position within +/-100 kb of the same TSS
            local = np.zeros((NSHUF, len(sub)), dtype=bool)
            for i, r in enumerate(sub.itertuples()):
                off = rng.integers(-100_000, 100_001, NSHUF)
                p = np.clip(int(r.tss) + off, 0, None)
                for j in range(NSHUF):
                    local[j, i] = overlap(iv, r.chrom_c, int(p[j]))
            local_frac = local.mean()

            rec = dict(cls=cls, n=len(sub), n_obs=int(obs.sum()),
                       observed_frac=obs.mean(),
                       local_null_frac=float(local_frac),
                       local_ratio=(obs.mean() / local_frac
                                    if local_frac > 0 else np.nan))
            k = obs.astype(np.int8) + local.sum(axis=0).astype(np.int8)
            kc = np.bincount(k, minlength=NSHUF + 2)
            rec["local_p"] = poisson_binomial_p(kc, int(obs.sum()), NSHUF)

            if cls in ("intronic", "antisense"):
                hobs, hnull = [], []
                for r in sub.itertuples():
                    h = r.partner_pcg
                    if not isinstance(h, str) or h not in hosts.index:
                        continue
                    hs, he = int(hosts.at[h, "start"]), int(hosts.at[h, "end"])
                    if he <= hs:
                        continue
                    hobs.append(overlap(iv, r.chrom_c, int(r.tss)))
                    picks = rng.integers(hs, he, NSHUF)
                    hnull.append([overlap(iv, hosts.at[h, "chrom_c"],
                                          int(p)) for p in picks])
                if len(hobs) >= 30:
                    hobs = np.array(hobs)
                    hnull = np.array(hnull, dtype=bool)
                    rec["host_n"] = len(hobs)
                    rec["host_observed_frac"] = hobs.mean()
                    rec["host_null_frac"] = hnull.mean()
                    rec["host_ratio"] = (hobs.mean() / hnull.mean()
                                         if hnull.mean() > 0 else np.nan)
                    # hnull is genes x shuffles, one row per host gene
                    k = hobs.astype(np.int8) + hnull.sum(axis=1).astype(np.int8)
                    kc = np.bincount(k, minlength=NSHUF + 2)
                    rec["host_p"] = poisson_binomial_p(kc, int(hobs.sum()),
                                                       NSHUF)
            rows.append(rec)

        T = pd.DataFrame(rows)
        # Direct class-versus-PCG contrast on the observed overlap counts
        from scipy.stats import fisher_exact
        pcg = T[T.cls == "PCG"]
        if len(pcg):
            pcg = pcg.iloc[0]
            T["vs_pcg_OR"] = np.nan
            T["vs_pcg_p"] = np.nan
            for i, r in T.iterrows():
                if r.cls == "PCG":
                    continue
                orr, p = fisher_exact([[int(r.n_obs), int(r.n - r.n_obs)],
                                       [int(pcg.n_obs),
                                        int(pcg.n - pcg.n_obs)]])
                T.loc[i, "vs_pcg_OR"] = orr
                T.loc[i, "vs_pcg_p"] = p
        for col, fdr in (("local_p", "local_fdr"), ("host_p", "host_fdr"),
                         ("vs_pcg_p", "vs_pcg_fdr")):
            if col in T:
                v = T[col].astype(float).values
                ok = ~np.isnan(v)
                f = np.full(len(v), np.nan)
                if ok.sum():
                    f[ok] = bh_fdr(v[ok])
                T[fdr] = f
        T.to_csv(safe_write_path(os.path.join(
            OUT, "E_erna_overlap.tsv")), sep="\t", index=False)

        # LTR contrast within intergenic, matched on eRNA neighbourhood
        ig = m[m.cls == "intergenic"].copy()
        ig["obs"] = [overlap(iv, r.chrom_c, int(r.tss)) for r in ig.itertuples()]
        a = ig[ig.tss_in_ltr == 1]; b = ig[ig.tss_in_ltr != 1]
        if len(a) > 20 and len(b) > 20:
            st.note("intergenic eRNA overlap LTR {:.1%} v other {:.1%}".format(
                a.obs.mean(), b.obs.mean()))

        for r in T.itertuples():
            st.note("{}: observed {:.1%} v local null {:.1%} (ratio {:.2f}, "
                    "P={:.2g}, n={})".format(
                        r.cls, r.observed_frac, r.local_null_frac,
                        r.local_ratio, r.local_p, r.n))
            if hasattr(r, "host_ratio") and not pd.isna(r.host_ratio):
                st.note("   host null: {:.1%} v {:.1%} (ratio {:.2f}, "
                        "P={:.2g}, n={})".format(
                            r.host_observed_frac, r.host_null_frac,
                            r.host_ratio, r.host_p, r.host_n))
        st.summary_md(os.path.join(
            RES, "step_summaries", "74_erna_overlap.md"),
            "74 class TSSs and eRNAbase eRNAs",
            df_to_md(T, ".4f"))


if __name__ == "__main__":
    main()
