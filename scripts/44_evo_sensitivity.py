#!/usr/bin/env python
"""44_evo_sensitivity.py — sensitivity analyses for the evolutionary layer.

1. Promoter-vs-body comparison with covariates (GC, gene length, expression,
   TE coverage): multi-covariate rank-residual partial Spearman, plus
   nearest-neighbour matched differences (class vs PCG).
2. High-confidence lncRNA subset (CAGE-refined TSS) sensitivity for key
   conservation/patch/TE metrics.
3. Age-assignment sensitivity: age-class composition and age-phyloP
   correlation in all vs high-confidence lncRNAs.
Outputs -> results/evo_sensitivity/*.tsv + step summary.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "evo_sensitivity")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
METRICS = ["prom_minus_body_phylop", "prom_minus_body_phast",
           "prom_phylop", "body_phylop", "exon_phylop",
           "prom_patch_frac_phylop", "te_exon_frac", "ltr_exon_frac"]
COVARS = ["gc", "log_body_len", "log_max_tpm", "te_cov"]


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


def partial_rank(x, y, cov):
    """Spearman partial correlation of x,y controlling multiple covariates."""
    d = pd.DataFrame({"x": x, "y": y})
    for i, c in enumerate(cov):
        d["z%d" % i] = c
    d = d.dropna()
    if len(d) < 30:
        return np.nan, np.nan, len(d)
    R = d.rank()
    X = R[["x"]].values
    Y = R[["y"]].values
    Z = np.c_[R[["z%d" % i for i in range(len(cov))]].values,
              np.ones(len(d))]
    bx, *_ = np.linalg.lstsq(Z, X, rcond=None)
    by, *_ = np.linalg.lstsq(Z, Y, rcond=None)
    ex, ey = (X - Z @ bx).ravel(), (Y - Z @ by).ravel()
    r, p = spearmanr(ex, ey)
    return float(r), float(p), len(d)


def main():
    with step("44_evo_sensitivity") as st:
        raw = pd.read_csv(os.path.join(RES, "constraint_structure",
                                       "gene_profiles_raw.tsv"),
                          sep="\t", index_col=0, low_memory=False)
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        seq = pd.read_csv(os.path.join(RES, "seq_features.tsv"), sep="\t",
                          index_col=0, low_memory=False)[["gc"]]
        act = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                          index_col=0, low_memory=False)[["max_tpm"]]
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                         index_col=0, low_memory=False)[["te_cov"]]
        cte = pd.read_csv(os.path.join(RES, "constraint_structure",
                                       "C_te_per_gene.tsv"), sep="\t",
                          index_col=0)[["te_exon_frac", "ltr_exon_frac"]]
        m = meta[["cls", "group", "cage_refined"]].join(raw).join(seq) \
            .join(act).join(te).join(cte)
        m["cls"] = np.where(m.group == "pcg", "PCG", m.cls)
        m["log_body_len"] = np.log1p(m.body_len.clip(lower=0))
        m["log_max_tpm"] = np.log1p(m.max_tpm.clip(lower=0))
        m["prom_minus_body_phylop"] = m.prom_phylop - m.body_phylop
        m["prom_minus_body_phast"] = m.prom_phast - m.body_phast

        # ---- 1a. partial correlations ---------------------------------
        rows = []
        pcg = m[m.group == "pcg"]
        for met in METRICS:
            for c in CLASSES:
                s = m[(m.group == "pcg") | (m.cls == c)]
                x = (s.cls == c).astype(int).values
                cov = [s[cv].values for cv in COVARS]
                r, p, n = partial_rank(x, s[met].values, cov)
                rows.append(dict(metric=met, cls=c, n=n,
                                 partial_rho=r, pval=p))
        P = pd.DataFrame(rows)
        P["fdr"] = bh_fdr(P.pval.fillna(1).values)
        P.to_csv(safe_write_path(os.path.join(
            OUT, "partial_corr.tsv")), sep="\t", index=False)

        # ---- 1b. nearest-neighbour matched differences ----------------
        rows = []
        for met in METRICS:
            pcg_met = pcg.dropna(subset=COVARS + [met])
            if len(pcg_met) < 100:
                continue
            scaler = StandardScaler().fit(pcg_met[COVARS])
            nn = NearestNeighbors(n_neighbors=5).fit(
                scaler.transform(pcg_met[COVARS]))
            for c in CLASSES:
                cls_df = m[(m.cls == c) & m[met].notna()]
                cls_df = cls_df.dropna(subset=COVARS)
                if len(cls_df) < 20:
                    continue
                dist, idx = nn.kneighbors(
                    scaler.transform(cls_df[COVARS]))
                ctrl = pcg_met[met].values[idx].mean(axis=1)
                diff = cls_df[met].values - ctrl
                rows.append(dict(metric=met, cls=c, n=len(cls_df),
                                 median_class=float(np.median(
                                     cls_df[met].values)),
                                 median_matched_PCG=float(np.median(ctrl)),
                                 median_diff=float(np.median(diff)),
                                 mean_diff=float(np.mean(diff))))
        M = pd.DataFrame(rows)
        M.to_csv(safe_write_path(os.path.join(
            OUT, "matched_diff.tsv")), sep="\t", index=False)

        # ---- 2. high-confidence lncRNA subset --------------------------
        hc = m[(m.group == "pcg") | (m.cage_refined == 1)]
        rows = []
        for met in METRICS:
            x = hc.loc[hc.group == "pcg", met].dropna().values
            for c in CLASSES:
                y = hc.loc[hc.cls == c, met].dropna().values
                if len(y) < 20:
                    continue
                d, p = mwu(y, x)
                rows.append(dict(metric=met, cls=c, n=len(y),
                                 median_class=float(np.median(y)),
                                 median_PCG=float(np.median(x)),
                                 cliffs_delta=d, pval=p))
        H = pd.DataFrame(rows)
        H["fdr"] = bh_fdr(H.pval.values)
        H.to_csv(safe_write_path(os.path.join(
            OUT, "highconf_tests.tsv")), sep="\t", index=False)

        # ---- 3. age assignment sensitivity ----------------------------
        age = pd.read_csv(os.path.join(RES, "comparative",
                                       "E_transcript_age.tsv"), sep="\t",
                          index_col=0)[["age_level", "age_class"]]
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)[
            ["phylop_mean"]]
        a = meta[["group", "cage_refined"]].join(age).join(ph)
        a = a[a.group == "lnc"]
        rows = []
        for label, sub in [("all_lnc", a),
                           ("highconf_lnc", a[a.cage_refined == 1])]:
            vc = sub.age_class.value_counts(normalize=True)
            ok = sub[["age_level", "phylop_mean"]].dropna()
            r, p = spearmanr(ok.age_level, ok.phylop_mean)
            rows.append(dict(group=label, n=len(sub),
                             human_specific=float(vc.get("human_specific", 0)),
                             primate=float(vc.get("primate", 0)),
                             mammalian=float(vc.get("mammalian", 0)),
                             eutherian=float(vc.get("eutherian", 0)),
                             rho_age_phylop=float(r), pval=float(p),
                             n_corr=len(ok)))
        A = pd.DataFrame(rows)
        A.to_csv(safe_write_path(os.path.join(
            OUT, "age_sensitivity.tsv")), sep="\t", index=False)

        st.note("partial corr: " + "; ".join(
            "{} {} rho={:.3f} p={:.1e}".format(r.metric, r.cls,
                                               r.partial_rho, r.pval)
            for r in P.itertuples() if r.metric ==
            "prom_minus_body_phylop"))
        st.note("matched diff: " + "; ".join(
            "{} {} med={:+.3f}".format(r.metric, r.cls, r.median_diff)
            for r in M.itertuples() if r.metric ==
            "prom_minus_body_phylop"))
        st.note("age sensitivity: " + "; ".join(
            "{} n={} human={:.1%} rho={:.3f}".format(
                r.group, r.n, r.human_specific, r.rho_age_phylop)
            for r in A.itertuples()))
        st.summary_md("results/evo_sensitivity/step_summaries/44.md",
                      "44 evolutionary sensitivity analyses",
                      "## Partial Spearman (GC, length, expression, TE)\n\n"
                      + df_to_md(P, ".3f", max_rows=60)
                      + "\n\n## Matched differences (class vs 5 PCG neighbours)\n\n"
                      + df_to_md(M, ".3f", max_rows=60)
                      + "\n\n## High-confidence lncRNA subset\n\n"
                      + df_to_md(H, ".3f", max_rows=60)
                      + "\n\n## Age sensitivity\n\n" + df_to_md(A, ".3f"))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    main()
