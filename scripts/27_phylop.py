#!/usr/bin/env python
"""27_phylop.py — phyloP 241-mammal conservation + acceleration layer.

Complements the phastCons470way layer (10_consv_chromatin) with per-base
phyloP from the Cactus 241-mammal alignment over the SAME -500/+200
promoter windows:
  A. class-vs-PCG tests on phylop_mean / frac_accel_2 / frac_cons_2
     (MWU + Cliff's delta + BH; own table — 12a pub_stats untouched)
  B. GC partial Spearman (12d rank-residual method) on phylop_mean
  C. acceleration / LTR-origin contrasts (frac_accel_2, phylop_mean)
  D. phastCons-vs-phyloP triangulation (class means, LTR subset,
     within-class cross-track Spearman)
Outputs: results/phylop/{phylop_per_gene,pub_stats_phylop,
     partial_corr_phylop,accel_tests_phylop,phastcons_vs_phylop}.tsv
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
PYLOP = os.path.join(RES, "phylop")
SIF_PBW = "/media/share/resource/container/pyBigWig/pyBigWig.sif"

CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic", "pseudogene", "transcribed_pseudogene"]
LNC5 = ["divergent", "antisense", "intronic", "sense_overlapping",
        "intergenic"]


def sing(*args):
    cmd = ["singularity", "exec", "--no-mount", "bind-paths",
           "--bind", "/media/share/resource:/media/share/resource:ro",
           "--bind", f"{WORKDIR}:{WORKDIR}"] + list(args)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(" ".join(cmd) + "\n" + r.stderr[:1000])
    return r.stdout


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def partial_spearman(x, y, z):
    """Spearman partial correlation of x,y controlling z (rank residuals)."""
    rx, ry, rz = pd.Series(x).rank(), pd.Series(y).rank(), pd.Series(z).rank()
    def resid(a):
        A = np.c_[rz, np.ones(len(rz))]
        b, *_ = np.linalg.lstsq(A, a.values, rcond=None)
        return a.values - A @ b
    ex, ey = resid(rx), resid(ry)
    r, p = spearmanr(ex, ey)
    return r, p


def mwu(y, x):
    """Two-sided MWU + Cliff's delta of y (class) vs x (reference)."""
    u, p = mannwhitneyu(y, x, alternative="two-sided")
    delta = 2 * u / (len(y) * len(x)) - 1
    return float(delta), float(p)


def main():
    with step("27_phylop") as st:
        # ---- container pass (per-base phyloP), skipped if present
        raw = os.path.join(PYLOP, "phylop_raw.tsv")
        if not os.path.exists(raw):
            sing(SIF_PBW, "python",
                 os.path.join(WORKDIR, "scripts", "pybw_stats_phylop.py"))
        ph = pd.read_csv(raw, sep="\t", index_col=0)
        st.note(f"phylop_raw rows={len(ph)}; nan mean="
                f"{int(ph.phylop_mean.isna().sum())}")

        # ---- minimal join (no 25-col master drift)
        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)[["group", "cls"]]
        m = m.join(ph, how="left")
        m = m.join(pd.read_csv(os.path.join(RES, "seq_features.tsv"),
                               sep="\t", index_col=0,
                               low_memory=False)[["gc"]])
        m = m.join(pd.read_csv(os.path.join(RES, "conservation.tsv"),
                               sep="\t", index_col=0,
                               low_memory=False)[["phast_mean"]])
        m = m.join(pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"),
                               sep="\t", index_col=0,
                               low_memory=False)[["tss_in_ltr"]])
        m.to_csv(safe_write_path(os.path.join(PYLOP, "phylop_per_gene.tsv")),
                 sep="\t")

        # ---- A. class vs PCG (12a machinery, own BH family)
        pcg = m[m.group == "pcg"]
        rows = []
        for feat, label in [("phylop_mean", "phyloP (241w) mean"),
                            ("frac_accel_2", "frac accelerated (phylop<=-2)"),
                            ("frac_cons_2", "frac conserved (phylop>=+2)")]:
            x = pcg[feat].dropna().values
            for cls in CLASSES:
                y = m.loc[m.cls == cls, feat].dropna().values
                if len(y) < 20:
                    continue
                delta, p = mwu(y, x)
                rows.append(dict(feature=label, cls=cls, n_class=len(y),
                                 mean_class=float(np.mean(y)),
                                 mean_PCG=float(np.mean(x)),
                                 cliffs_delta=delta, pval=p))
        A = pd.DataFrame(rows)
        A["fdr"] = bh_fdr(A["pval"].values)
        A = A.sort_values(["feature", "cls"])
        A.to_csv(safe_write_path(
            os.path.join(PYLOP, "pub_stats_phylop.tsv")), sep="\t",
            index=False)
        st.note(f"A: rows={len(A)}; significant(FDR<0.05)="
                f"{int((A.fdr < 0.05).sum())}")

        # ---- B. GC partial correlation (12d section-3 mirror)
        prows = []
        sub = m[m.group.isin(["pcg", "lnc"]) & m.phylop_mean.notna()
                & m.gc.notna()]
        for cls in LNC5:
            s = sub[(sub.group == "pcg") | (sub.cls == cls)].copy()
            x = (s["cls"] == cls).astype(int).values
            r0, p0 = spearmanr(x, s["phylop_mean"].values)
            r, p = partial_spearman(pd.Series(x), s["phylop_mean"], s["gc"])
            prows.append(dict(cls=cls, n=len(s), rho_raw=r0, p_raw=p0,
                              rho_GC_controlled=r, p_controlled=p))
        B = pd.DataFrame(prows)
        B.to_csv(safe_write_path(
            os.path.join(PYLOP, "partial_corr_phylop.tsv")), sep="\t",
            index=False)
        st.note("B phyloP GC-controlled: " + "; ".join(
            f"{r.cls}: {r.rho_raw:.3f}->{r.rho_GC_controlled:.3f}"
            for r in B.itertuples()))

        # ---- C. acceleration / LTR-origin contrasts
        crows = []
        contrasts = [
            ("lnc LTR-TSS vs non-LTR",
             (m.group == "lnc") & (m.tss_in_ltr == 1),
             (m.group == "lnc") & (m.tss_in_ltr == 0)),
            ("intergenic LTR-TSS vs non-LTR",
             (m.cls == "intergenic") & (m.tss_in_ltr == 1),
             (m.cls == "intergenic") & (m.tss_in_ltr == 0))]
        for name, cmask, kmask in contrasts:
            for feat in ["frac_accel_2", "phylop_mean"]:
                y = m.loc[cmask, feat].dropna().values
                x = m.loc[kmask, feat].dropna().values
                if len(y) < 20 or len(x) < 20:
                    continue
                delta, p = mwu(y, x)
                crows.append(dict(contrast=name, feature=feat,
                                  n_case=len(y), n_ctrl=len(x),
                                  median_case=float(np.median(y)),
                                  median_ctrl=float(np.median(x)),
                                  mean_case=float(np.mean(y)),
                                  mean_ctrl=float(np.mean(x)),
                                  cliffs_delta=delta, pval=p))
        C = pd.DataFrame(crows)
        C["fdr"] = bh_fdr(C["pval"].values)
        C.to_csv(safe_write_path(
            os.path.join(PYLOP, "accel_tests_phylop.tsv")), sep="\t",
            index=False)
        st.note("C LTR contrasts: " + "; ".join(
            f"{r.contrast.split()[0]}/{r.feature}={r.mean_case:.4f} v "
            f"{r.mean_ctrl:.4f} (p={r.pval:.2e})"
            for r in C[C.feature == "frac_accel_2"].itertuples()))

        # ---- D. phastCons vs phyloP triangulation
        drows = []
        groups = [("PCG", m.group == "pcg")] + [(c, m.cls == c)
                                                for c in CLASSES]
        for name, mask in groups + [("all_promoters", m.index == m.index)]:
            s = m[mask]
            ok = s[["phast_mean", "phylop_mean"]].dropna()
            if len(ok) >= 20:
                r, p = spearmanr(ok["phast_mean"], ok["phylop_mean"])
            else:
                r, p = np.nan, np.nan
            ltr = s[s.tss_in_ltr == 1]
            drows.append(dict(
                cls=name, n=len(s),
                phast_mean=float(s.phast_mean.mean()),
                phylop_mean=float(s.phylop_mean.mean()),
                frac_accel_2=float(s.frac_accel_2.mean()),
                frac_cons_2=float(s.frac_cons_2.mean()),
                n_ltr=int(len(ltr)),
                frac_accel2_ltr=float(ltr.frac_accel_2.mean())
                if len(ltr) >= 20 else float("nan"),
                spearman_phast_phylop=float(r), p_spearman=float(p)))
        D = pd.DataFrame(drows)
        D.to_csv(safe_write_path(
            os.path.join(PYLOP, "phastcons_vs_phylop.tsv")), sep="\t",
            index=False)
        st.note("D cross-track spearman(all)="
                f"{D.loc[D.cls == 'all_promoters', 'spearman_phast_phylop'].iloc[0]:.3f}")

        body = ("## A. Class vs PCG (phyloP metrics)\n\n"
                + df_to_md(A, ".4f", max_rows=25)
                + "\n\n## B. GC-controlled partial Spearman (phylop_mean)\n\n"
                + df_to_md(B, ".3f")
                + "\n\n## C. LTR-origin acceleration contrasts\n\n"
                + df_to_md(C, ".4f")
                + "\n\n## D. phastCons vs phyloP by class\n\n"
                + df_to_md(D, ".3f"))
        st.summary_md("results/step_summaries/27_phylop.md",
                      "27 phyloP 241-mammal conservation + acceleration",
                      body)


if __name__ == "__main__":
    main()
