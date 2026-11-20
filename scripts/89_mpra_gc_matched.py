#!/usr/bin/env python
"""89_mpra_gc_matched.py — GC-matched re-analysis of the ENCODE lentiMPRA.

The element-level comparison in results/mpra compares LTR-overlapping DNase
elements with all other elements. LTRs are GC-rich and GC drives reporter
activity, so the difference is an upper bound. This module recomputes the
contrast within GC bins using the local hg38 FASTA, and repeats it for elements
that overlap an LTR-donated lncRNA TSS. The assay still uses a minimal
promoter, so the readout is regulatory potential rather than promoter autonomy.

Outputs
  results/mpra/G_gc_bins.tsv
  results/mpra/G_summary.tsv
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
import pysam
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "mpra")
ELEM = os.path.join(OUT, "raw_encode_lentimpra", "ENCFF613VFY.bed.gz")
LTR = os.path.join(WORKDIR, "data", "constraint_structure_tmp", "TE_LTR.bed")
FA = "/media/share/resource/reference/homo_sapiens/GRCh38/genome/GRCh38.CHR.fa"
REF = os.path.join(OUT, "ltr_tss_reference.bed")


def load_elements():
    rows = []
    with gzip.open(ELEM, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            try:
                act = float(f[6])
            except ValueError:
                continue
            rows.append((f[0], int(f[1]), int(f[2]), f[3], act))
    d = pd.DataFrame(rows, columns=["chrom", "start", "end", "name", "act"])
    g = d.groupby(["chrom", "start", "end"], as_index=False).agg(
        act=("act", "mean"), name=("name", "first"))
    return g


def load_intervals(path, cols):
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < cols[1] + 1:
                continue
            try:
                a, b = int(f[cols[0]]), int(f[cols[1]])
            except ValueError:
                continue
            out.setdefault(f[0], []).append((a, b))
    for c in out:
        out[c] = np.array(sorted(out[c]))
    return out


def overlap_frac(d, iv):
    """Fraction of each element covered by intervals in iv."""
    frac = np.zeros(len(d))
    for c, g in d.groupby("chrom"):
        t = iv.get(c)
        if t is None:
            continue
        s, e = t[:, 0], t[:, 1]
        for i, r in zip(g.index, g.itertuples()):
            j0 = np.searchsorted(e, r.start, "right")
            j1 = np.searchsorted(s, r.end, "left")
            if j1 <= j0:
                continue
            lo = max(r.start, s[j0])
            hi = min(r.end, e[j1 - 1])
            cov = max(0, hi - lo)
            frac[i] = min(1.0, cov / float(r.end - r.start))
    return frac


def main():
    with step("89_mpra_gc_matched") as st:
        fa = pysam.FastaFile(FA)
        d = load_elements()
        st.note("unique elements={}".format(len(d)))
        gc = np.zeros(len(d))
        for i, r in enumerate(d.itertuples()):
            try:
                s = fa.fetch(r.chrom, r.start, r.end).upper()
            except (KeyError, ValueError):
                continue
            gc[i] = (s.count("G") + s.count("C")) / max(1, len(s))
        d["gc"] = gc
        d = d[d.gc > 0].copy().reset_index(drop=True)

        ltr_iv = load_intervals(LTR, (1, 2))
        d["ltr_frac"] = overlap_frac(d, ltr_iv)
        tss = pd.read_csv(REF, sep="\t", header=None, comment="#")
        tss_iv = {}
        for r in tss.itertuples():
            tss_iv.setdefault(r[1], []).append((int(r[5]) - 100,
                                                int(r[5]) + 100))
        for c in tss_iv:
            tss_iv[c] = np.array(sorted(tss_iv[c]))
        d["tss_frac"] = overlap_frac(d, tss_iv)

        dnase = d[d.name.str.startswith("DNasePeak")].copy()
        dnase["gcbin"] = pd.cut(dnase.gc, np.arange(0.30, 0.86, 0.02))
        dnase["is_ltr"] = dnase.ltr_frac > 0

        def strat(case, ctrl, label):
            rows = []
            for b in dnase.gcbin.cat.categories:
                a = case.act[case.gcbin == b]
                n = ctrl.act[ctrl.gcbin == b]
                if len(a) < 20 or len(n) < 20:
                    continue
                rows.append(dict(subset=label, gc_bin=str(b),
                                 n_ltr=len(a), n_non=len(n),
                                 mean_ltr=float(a.mean()),
                                 mean_non=float(n.mean()),
                                 delta=float(a.mean() - n.mean()),
                                 frac_active_ltr=float((a > 0).mean()),
                                 frac_active_non=float((n > 0).mean())))
            return pd.DataFrame(rows)

        B = strat(dnase[dnase.is_ltr], dnase[~dnase.is_ltr], "dnase_all")
        strict = dnase[dnase.ltr_frac >= 0.5]
        B2 = strat(strict, dnase[~dnase.is_ltr], "ltr_half_or_more")
        B = pd.concat([B, B2], ignore_index=True)
        B.to_csv(safe_write_path(os.path.join(OUT, "G_gc_bins.tsv")),
                 sep="\t", index=False)

        def wdelta(t):
            w = (t.n_ltr + t.n_non).values
            return float((t.delta.values * w).sum() / w.sum())

        delta = wdelta(B[B.subset == "dnase_all"])
        delta_strict = wdelta(B2) if len(B2) else np.nan

        # permutation: shuffle the LTR label within GC bin, 2000 draws
        rng = np.random.default_rng(42)
        sub = dnase[dnase.gcbin.isin(B[B.subset == "dnase_all"].gc_bin.astype(
            str).values)].copy()
        codes = sub.gcbin.cat.codes.values
        lab = sub.is_ltr.values.astype(bool)
        act = sub.act.values
        wts = {}
        for c in np.unique(codes):
            wts[c] = (sub.gcbin.cat.categories[c] in
                      set(B[B.subset == "dnase_all"].gc_bin))
        cnt = 0
        nperm = 2000
        for _ in range(nperm):
            pl = lab.copy()
            for c in np.unique(codes):
                i = np.where(codes == c)[0]
                pl[i] = rng.permutation(pl[i])
            num = den = 0.0
            for c in np.unique(codes):
                i = np.where(codes == c)[0]
                a = act[i][pl[i]]
                n = act[i][~pl[i]]
                if len(a) < 20 or len(n) < 20:
                    continue
                w = len(a) * len(n) / (len(a) + len(n))
                num += (a.mean() - n.mean()) * w
                den += w
            if den and abs(num / den) >= abs(delta):
                cnt += 1
        p_perm = (cnt + 1) / (nperm + 1)

        S = pd.DataFrame([
            dict(subset="dnase_all", n_ltr=int(dnase.is_ltr.sum()),
                 n_non=int((~dnase.is_ltr).sum()),
                 raw_delta=float(dnase[dnase.is_ltr].act.mean()
                                 - dnase[~dnase.is_ltr].act.mean()),
                 gc_matched_delta=delta, perm_p=p_perm),
            dict(subset="ltr_half_or_more", n_ltr=len(strict),
                 n_non=int((~dnase.is_ltr).sum()),
                 raw_delta=float(strict.act.mean()
                                 - dnase[~dnase.is_ltr].act.mean()),
                 gc_matched_delta=delta_strict, perm_p=np.nan),
            dict(subset="ltr_tss_100bp",
                 n_ltr=int((dnase.tss_frac > 0).sum()),
                 n_non=int((dnase.tss_frac == 0).sum()),
                 raw_delta=float(dnase[dnase.tss_frac > 0].act.mean()
                                 - dnase[dnase.tss_frac == 0].act.mean()),
                 gc_matched_delta=np.nan, perm_p=np.nan)])
        S.to_csv(safe_write_path(os.path.join(OUT, "G_summary.tsv")),
                 sep="\t", index=False)

        st.note("GC range LTR {:.3f} v non-LTR {:.3f}".format(
            dnase[dnase.ltr_frac > 0].gc.mean(),
            dnase[dnase.ltr_frac == 0].gc.mean()))
        st.note("raw delta {:.3f}, GC-matched delta {:.3f}, perm P={:.3f}"
                .format(S.raw_delta.iloc[0], delta, p_perm))
        st.note("elements at LTR TSS (+/-100 bp): {} of {}".format(
            int((dnase.tss_frac > 0).sum()), len(dnase)))
        st.summary_md(os.path.join(
            RES, "step_summaries", "89_mpra_gc_matched.md"),
            "89 GC-matched lentiMPRA activity",
            df_to_md(B, ".3f") + "\n\n" + df_to_md(S, ".3f"))


if __name__ == "__main__":
    main()
