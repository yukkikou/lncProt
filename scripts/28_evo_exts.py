#!/usr/bin/env python
"""28_evo_exts.py — three evolution-extension analyses (analysis-first round;
manuscript untouched, results only).

Modes (13d-style CLI, each restartable via its own step marker):
  A -> ML feature ablation: does phyloP add classification accuracy for
       birth class beyond phastCons? Feature sets F0 (base-19, must reproduce
       module 24 macro-F1 0.336), F1 (-phast), F2 (+5 phyloP), F3 (phyloP
       replaces phast). Same seeded folds; seeded permutation importance.
  B -> broad vs tissue-specific lncRNAs: class composition + evolution
       metrics (GTEx tau quartiles primary, CAGE breadth replication).
  C -> cancer-associated lncRNAs: per-gene TCGA DE recomputed (13d logic,
       bit-compatible spot-check vs shards), up-tier vs class/evolution.
  R -> append three rows to REPORT.md appendix.
Outputs -> results/evo_ext/*.tsv + results/evo_ext/step_summaries/28_*.md
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu, spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import f1_score
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "evo_ext")
TD = os.path.join(WORKDIR, "data", "tissue_disease")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
PHYLOP5 = ["phylop_mean", "phylop_max", "phylop_min", "frac_cons_2",
           "frac_accel_2"]
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


def mwu_row(y, x):
    u, p = mannwhitneyu(y, x, alternative="two-sided")
    return float(p), float(np.median(y)), float(np.median(x))


# ------------------------------------------------------------------ mode A --
def build_feature_table():
    """Verbatim from 24_ml_extensions.build_feature_table (L39-61)."""
    m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                    index_col=0, low_memory=False)
    parts = [m[["cls", "group"]]]
    for f, cols in [("seq_features.tsv", ["gc", "cpg_oe", "cpg_density",
                                          "TATA", "Inr", "CCAAT", "GCbox",
                                          "cpg_island", "prom_class"]),
                    ("te_promoter_stats.tsv", ["te_cov", "te_frac_LTR",
                                               "te_frac_LINE", "te_frac_SINE",
                                               "te_frac_DNA", "tss_in_te",
                                               "tss_in_ltr"]),
                    ("conservation.tsv", ["phast_mean"])]:
        d = pd.read_csv(os.path.join(RES, f), sep="\t", index_col=0,
                        low_memory=False)
        parts.append(d[[c for c in cols if c in d.columns]])
    ch = pd.read_csv(os.path.join(RES, "chromhmm_states.tsv"), sep="\t",
                     index_col=0, low_memory=False)
    prom_like = ch.eq("TssA").mean(axis=1).rename("prom_like")
    enh_like = ch.isin(["EnhG", "Enh"]).mean(axis=1).rename("enh_like")
    parts.append(pd.concat([prom_like, enh_like], axis=1))
    X = pd.concat(parts, axis=1)
    X["prom_class"] = X["prom_class"].map({"HCP": 2, "ICP": 1, "LCP": 0})
    return X


def run_cv(Xs, y):
    """Verbatim CV block from 24 task_a (no confusion output)."""
    skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
    clf = HistGradientBoostingClassifier(max_depth=6, learning_rate=0.08,
                                         max_iter=400, random_state=SEED)
    yhat = np.empty(len(y), dtype=object)
    for tr, te in skf.split(Xs, y):
        clf.fit(Xs.iloc[tr], y[tr])
        yhat[te] = clf.predict(Xs.iloc[te])
    return f1_score(y, yhat, average="macro"), float((y == yhat).mean())


def mode_a():
    with step("28_evo_A") as st:
        X = build_feature_table()
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)
        X = X.join(ph[PHYLOP5])
        sub = X[X.cls.isin(CLASSES)].copy()
        y = sub.cls.values
        base_cols = [c for c in sub.columns if c not in
                     ("cls", "group") + tuple(PHYLOP5)]
        base_feats = [c for c in base_cols if sub[c].notna().mean() > 0.8]
        phy_feats = [c for c in PHYLOP5 if sub[c].notna().mean() > 0.8]

        sets = [("F0_base19", base_feats),
                ("F1_minus_phast", [c for c in base_feats
                                    if c != "phast_mean"]),
                ("F2_plus_phylop", base_feats + phy_feats),
                ("F3_phylop_replaces_phast",
                 [c for c in base_feats if c != "phast_mean"] + phy_feats)]
        rows = []
        f0_macro = None
        for name, feats in sets:
            Xs = sub[feats].fillna(sub[feats].median(numeric_only=True))
            macro, acc = run_cv(Xs, y)
            if name == "F0_base19":
                f0_macro = macro
                if round(macro, 3) != 0.336:      # hard determinism gate
                    raise RuntimeError(
                        "F0 gate FAILED: macro-F1 {} != 0.336 "
                        "(module 24 baseline)".format(macro))
            rows.append(dict(feature_set=name, n_features=len(feats),
                             macro_f1=macro, acc=acc,
                             delta_vs_F0=macro - f0_macro if f0_macro
                             else np.nan))
        A = pd.DataFrame(rows)
        A.to_csv(safe_write_path(os.path.join(OUT, "A_ablation.tsv")),
                 sep="\t", index=False)
        st.note("A ablation: " + "; ".join(
            "{}: F1={:.4f} (d={:+.4f})".format(r.feature_set, r.macro_f1,
                                               r.delta_vs_F0)
            for r in A.itertuples()))

        # seeded permutation importance on F2 (single holdout refit, 24 design)
        feats = base_feats + phy_feats
        Xs = sub[feats].fillna(sub[feats].median(numeric_only=True))
        tr, te = next(StratifiedKFold(5, shuffle=True,
                                      random_state=SEED).split(Xs, y))
        clf = HistGradientBoostingClassifier(max_depth=6, learning_rate=0.08,
                                             max_iter=400, random_state=SEED)
        clf.fit(Xs.iloc[tr], y[tr])
        base = f1_score(y[te], clf.predict(Xs.iloc[te]), average="macro")
        rng = np.random.default_rng(SEED)
        imps = []
        for c in feats:
            Xp = Xs.iloc[te].copy()
            Xp[c] = rng.permutation(Xp[c].values)
            imps.append(dict(feature=c, delta_f1=base - f1_score(
                y[te], clf.predict(Xp), average="macro")))
        imp = pd.DataFrame(imps).sort_values("delta_f1", ascending=False)
        imp.to_csv(safe_write_path(os.path.join(
            OUT, "A_feature_importance_F2.tsv")), sep="\t", index=False)
        st.note("F0 gate passed ({}); F2 top features: ".format(f0_macro)
                + ", ".join("{}({:+.3f})".format(r.feature, r.delta_f1)
                            for r in imp.head(6).itertuples()))

        st.summary_md("results/evo_ext/step_summaries/28_a.md",
                      "28A ML feature ablation",
                      "## Ablation (same 5-fold seed-42 folds)\n\n"
                      + df_to_md(A, ".4f")
                      + "\n\n## F2 seeded permutation importance\n\n"
                      + df_to_md(imp, ".4f", max_rows=30))


# ------------------------------------------------------------------ mode B --
def metric_battery(df, scope, rows):
    """MWU on continuous metrics + Fisher on tss_in_ltr, broad vs specific."""
    b = df[df._grp == "broad"]
    s = df[df._grp == "specific"]
    for m in ["phylop_mean", "phast_mean", "gc", "te_cov"]:
        y, x = b[m].dropna(), s[m].dropna()
        if len(y) >= 20 and len(x) >= 20:
            p, mb, ms = mwu_row(y.values, x.values)
            rows.append(dict(analysis="", scope=scope, metric=m, test="MWU",
                             pval=p, n_broad=len(y), n_specific=len(x),
                             median_broad=mb, median_specific=ms))
    a = int((b.tss_in_ltr == 1).sum()); c = int(len(b) - a)
    d = int((s.tss_in_ltr == 1).sum()); e = int(len(s) - d)
    orr, p = fisher_exact([[a, c], [d, e]])
    rows.append(dict(analysis="", scope=scope, metric="tss_in_ltr",
                     test="Fisher_OR", pval=p, n_broad=len(b),
                     n_specific=len(s), median_broad=a / max(len(b), 1),
                     median_specific=d / max(len(s), 1)))


def split_quartiles(df, col, lo, hi):
    df = df.copy()
    df["_grp"] = np.where(df[col] <= lo, "broad",
                          np.where(df[col] >= hi, "specific", "mid"))
    return df[df._grp != "mid"]


def mode_b():
    with step("28_evo_B") as st:
        gt = pd.read_csv(os.path.join(RES, "tissue_disease", "tables",
                                      "gtex_gene_tau.tsv"), sep="\t",
                         index_col=0)
        tg = pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t",
                         index_col=0)
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                         index_col=0, low_memory=False)
        ver = gt.join(tg[["gene_id"]])
        B = ver.set_index("gene_id").join(
            ph[["phylop_mean", "phast_mean", "gc", "tss_in_ltr"]]
        ).join(te[["te_cov"]])

        rows = []
        comp_rows = []
        # GTEx tau: lnc primary, PCG comparator (own quantiles)
        lnc_g = B[B.cls.isin(CLASSES) & B.tau_gtex.notna()]
        q25, q75 = lnc_g.tau_gtex.quantile([.25, .75])
        st.note("B gtex lnc n={} q25={:.3f} q75={:.3f}".format(
            len(lnc_g), q25, q75))
        lo, hi = lnc_g.tau_gtex.quantile([.25, .75])
        df = split_quartiles(lnc_g, "tau_gtex", lo, hi)
        tab = pd.crosstab(df.cls, df._grp).reindex(CLASSES, fill_value=0)
        chi2, pchi, _, _ = chi2_contingency(tab.values)
        rows.append(dict(analysis="gtex_lnc", scope="composition",
                         metric="cls x group", test="chi2", pval=pchi,
                         n_broad=int(tab.broad.sum()),
                         n_specific=int(tab.specific.sum()),
                         median_broad=np.nan, median_specific=np.nan))
        for c in CLASSES:
            comp_rows.append(dict(analysis="gtex_lnc", cls=c,
                                  n_broad=int(tab.loc[c, "broad"]),
                                  n_specific=int(tab.loc[c, "specific"]),
                                  frac_broad=tab.loc[c, "broad"]
                                  / max(int(tab.loc[c].sum()), 1)))
        metric_battery(df, "gtex_lnc_overall", rows)
        for c in CLASSES:
            sub = df[df.cls == c]
            if len(sub[sub._grp == "broad"]) >= 20 and \
                    len(sub[sub._grp == "specific"]) >= 20:
                metric_battery(sub, "gtex_lnc_" + c, rows)
        pcg_g = B[(B.cls == "PCG") & B.tau_gtex.notna()]
        lo, hi = pcg_g.tau_gtex.quantile([.25, .75])
        dfp = split_quartiles(pcg_g, "tau_gtex", lo, hi)
        metric_battery(dfp, "gtex_pcg", rows)
        # CAGE breadth replication (versioned ids, direct join)
        cage = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                           index_col=0, low_memory=False)
        cg = cage.join(ph[["phylop_mean", "phast_mean", "gc", "tss_in_ltr"]],
                       how="left").join(te[["te_cov"]], how="left")
        cg = cg[cg.cls.isin(CLASSES) & cg.breadth.notna()]
        lo, hi = cg.breadth.quantile([.25, .75])
        df = split_quartiles(cg, "breadth", lo, hi)
        st.note("B cage lnc n={} q25={:.0f} q75={:.0f}".format(
            len(cg), lo, hi))
        tab = pd.crosstab(df.cls, df._grp).reindex(CLASSES, fill_value=0)
        chi2, pchi, _, _ = chi2_contingency(tab.values)
        rows.append(dict(analysis="cage_lnc", scope="composition",
                         metric="cls x group", test="chi2", pval=pchi,
                         n_broad=int(tab.broad.sum()),
                         n_specific=int(tab.specific.sum()),
                         median_broad=np.nan, median_specific=np.nan))
        metric_battery(df, "cage_lnc_overall", rows)
        for c in CLASSES:
            sub = df[df.cls == c]
            if len(sub[sub._grp == "broad"]) >= 20 and \
                    len(sub[sub._grp == "specific"]) >= 20:
                metric_battery(sub, "cage_lnc_" + c, rows)

        T = pd.DataFrame(rows)
        T["fdr_bh"] = bh_fdr(T.pval.values)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "B_broad_vs_specific.tsv")), sep="\t", index=False)
        Cmp = pd.DataFrame(comp_rows)
        Cmp.to_csv(safe_write_path(os.path.join(OUT, "B_by_class.tsv")),
                   sep="\t", index=False)
        key = T[(T.scope.str.endswith("overall")) &
                (T.metric.isin(["phylop_mean", "tss_in_ltr"]))]
        st.note("B key tests: " + "; ".join(
            "{} {}: p={:.2e}".format(r.scope, r.metric, r.pval)
            for r in key.itertuples()))
        st.summary_md("results/evo_ext/step_summaries/28_b.md",
                      "28B broad vs tissue-specific lncRNAs",
                      "## Tests (BH within family)\n\n"
                      + df_to_md(T, ".4g", max_rows=60)
                      + "\n\n## GTEx class composition (broad vs specific)\n\n"
                      + df_to_md(Cmp, ".3f")
                      + "\n\nNote: tau split is WITHIN GTEx-expressing "
                        "lncRNAs (already globally specific, q25=0.834); "
                      "sense_overlapping cells exploratory (n=110).")


# ------------------------------------------------------------------ mode C --
def mode_c():
    with step("28_evo_C") as st:
        d = np.nan_to_num(np.load(os.path.join(TD, "toil_class_expr.npz"))
                          ["M"])
        genes = pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t",
                            index_col=0)
        samp = pd.read_csv(os.path.join(TD, "toil_samples.tsv"), sep="\t")
        lnc_mask = genes.cls.isin(CLASSES).values
        n_genes = len(genes)
        n_tested = np.zeros(n_genes, int)
        n_up = np.zeros(n_genes, int)
        n_down = np.zeros(n_genes, int)
        max_delta_up = np.full(n_genes, np.nan)
        min_fdr = np.full(n_genes, np.nan)

        cancers = []
        for cancer, sub in samp[samp.cancer.notna()].groupby("cancer"):
            if ((sub.condition == "tumor").sum() >= 20 and
                    (sub.condition == "normal").sum() >= 5):
                cancers.append(cancer)
        st.note("C cancers with >=20 tumor & >=5 normal: {} ({})".format(
            len(cancers), ",".join(c.replace("TCGA-", "") for c in cancers)))

        for cancer in cancers:
            tum = samp[(samp.cancer == cancer) &
                       (samp.condition == "tumor")].index.values
            nor = samp[(samp.cancer == cancer) &
                       (samp.condition == "normal")].index.values
            T = d[:, tum].mean(axis=1)
            N = d[:, nor].mean(axis=1)
            detect = (T > 0) | (N > 0)
            delta = T - N
            delta[~detect] = np.nan
            pvals = np.full(n_genes, np.nan)
            for i in range(n_genes):
                a = d[i, tum]
                b = d[i, nor]
                if a.std() == 0 and b.std() == 0:
                    continue
                _, p = mannwhitneyu(a, b, alternative="two-sided")
                pvals[i] = p
            fdr = bh_fdr(pvals[np.isfinite(pvals)])
            fdr_full = np.full(n_genes, np.nan)
            fdr_full[np.isfinite(pvals)] = fdr
            up = (delta > 0.5) & (fdr_full < 0.05)
            dn = (delta < -0.5) & (fdr_full < 0.05)

            # hard gate: class aggregates must match the 13d shard bit-wise
            shard = pd.read_csv(os.path.join(TD, "shards",
                                             "13d_{}.tsv".format(cancer)),
                                sep="\t").set_index("cls")
            for c in CLASSES + ["PCG"]:
                idx = (genes.cls == "PCG").values if c == "PCG" else \
                    (genes.cls == c).values
                den = int((idx & detect).sum())
                fu = float((up & idx).sum() / max(den, 1))
                fd = float((dn & idx).sum() / max(den, 1))
                med = float(np.nanmedian(delta[idx]))
                r = shard.loc[c]
                if (abs(fu - r.frac_up) > 1e-9 or
                        abs(fd - r.frac_down) > 1e-9 or
                        abs(med - r.median_delta) > 1e-9):
                    raise RuntimeError(
                        "13d spot-check FAILED for {} {} ({} vs shard "
                        "{}/{}/{})".format(cancer, c, fu, fd, med,
                                           r.frac_up, r.frac_down,
                                           r.median_delta))
            # accumulate per-gene (all genes; lnc subset written out)
            tested = np.isfinite(pvals)
            n_tested += (tested & detect).astype(int)
            n_up += (up).astype(int)
            n_down += (dn).astype(int)
            cand_max = np.where(up, delta, -np.inf)
            m = np.fmax(max_delta_up, cand_max)
            max_delta_up = np.where(np.isfinite(m), m, np.nan)
            min_fdr = np.fmin(min_fdr,
                              np.where(tested, fdr_full, np.nan))
            st.note("{}: tumor={} normal={} up(lnc)={}".format(
                cancer, len(tum), len(nor), int((up & lnc_mask).sum())))

        tier = np.where(n_up >= 5, 3,
                        np.where(n_up >= 2, 2, np.where(n_up == 1, 1, 0)))
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)
        # index by VERSIONED gene_id so the phylop join lands
        out = pd.DataFrame(index=pd.Index(genes["gene_id"].values,
                                          name="gene_id"))
        out["cls"] = genes.cls.values
        out["n_tested"] = n_tested
        out["n_up"] = n_up
        out["n_down"] = n_down
        out["max_delta_up"] = max_delta_up
        out["min_fdr"] = min_fdr
        out["tier"] = tier
        out = out.join(ph[["phylop_mean", "phast_mean", "gc", "tss_in_ltr"]])
        lnc_out = out[lnc_mask]
        lnc_out.to_csv(safe_write_path(os.path.join(
            OUT, "C_pergene_de_lnc.tsv")), sep="\t")

        rows = []
        for label, frame in [("", lnc_out),
                             ("_tested10", lnc_out[lnc_out.n_tested >= 10])]:
            tab = pd.crosstab(frame.tier, frame.cls)
            chi2, p, _, _ = chi2_contingency(tab.values)
            rows.append(dict(scope="composition" + label, metric="tier x cls",
                             test="chi2", pval=p,
                             n_group=int(tab.values.sum())))
            t3 = frame[frame.tier == 3]
            t0 = frame[frame.tier == 0]
            for m in ["phylop_mean", "phast_mean"]:
                if len(t3) >= 20 and len(t0) >= 20:
                    pv, m3, m0 = mwu_row(t3[m].dropna().values,
                                         t0[m].dropna().values)
                    rows.append(dict(scope="tier3_v_tier0" + label,
                                     metric=m, test="MWU", pval=pv,
                                     n_group=len(t3) + len(t0)))
                rho, prho = spearmanr(frame.tier, frame[m], nan_policy="omit")
                rows.append(dict(scope="spearman_tier" + label, metric=m,
                                 test="spearman", pval=prho, n_group=len(frame)))
            a = int((t3.tss_in_ltr == 1).sum()); c = int(len(t3) - a)
            d_ = int((t0.tss_in_ltr == 1).sum()); e = int(len(t0) - d_)
            orr, pf = fisher_exact([[a, c], [d_, e]])
            rows.append(dict(scope="tier3_v_tier0" + label,
                             metric="tss_in_ltr_frac", test="Fisher_OR",
                             pval=pf, n_group=len(t3) + len(t0)))
        # intergenic LTR x tier
        ig = lnc_out[lnc_out.cls == "intergenic"].copy()
        ig["tss_in_ltr"] = ig.tss_in_ltr.fillna(0).astype(int)
        tab = pd.crosstab(ig.tier, ig.tss_in_ltr)
        chi2, p, _, _ = chi2_contingency(tab.values)
        rows.append(dict(scope="intergenic_LTR_x_tier",
                         metric="tier x tss_in_ltr", test="chi2", pval=p,
                         n_group=len(ig)))
        T = pd.DataFrame(rows)
        T["fdr_bh"] = bh_fdr(T.pval.values)
        T.to_csv(safe_write_path(os.path.join(OUT, "C_tier_tests.tsv")),
                 sep="\t", index=False)
        tc = lnc_out.tier.value_counts().sort_index()
        st.note("C tiers (lnc n={}): ".format(len(lnc_out))
                + "; ".join("tier{}={}".format(t, int(n))
                            for t, n in tc.items())
                + "; spot-check 18x6 vs 13d shards PASSED")
        st.summary_md("results/evo_ext/step_summaries/28_c.md",
                      "28C cancer-associated lncRNAs (TCGA per-gene DE)",
                      "## Tier tests (BH within family)\n\n"
                      + df_to_md(T, ".4g")
                      + "\n\n## Tier x class counts\n\n"
                      + df_to_md(pd.crosstab(lnc_out.tier, lnc_out.cls), "d")
                      + "\n\nTier: 0=not_up_any, 1=up in 1, 2=up in 2-4, "
                        "3=up in >=5 of 18 cancers (n_tested>=10 rows are "
                        "the detection-confound sensitivity).")


# ------------------------------------------------------------------ mode R --
def mode_r():
    with step("28_evo_R") as st:
        A = pd.read_csv(os.path.join(OUT, "A_ablation.tsv"), sep="\t")
        f2 = A[A.feature_set == "F2_plus_phylop"].iloc[0]
        f3 = A[A.feature_set == "F3_phylop_replaces_phast"].iloc[0]
        B = pd.read_csv(os.path.join(OUT, "B_broad_vs_specific.tsv"), sep="\t")
        bl = B[(B.scope == "gtex_lnc_overall") &
               (B.metric == "phylop_mean")].iloc[0]
        C = pd.read_csv(os.path.join(OUT, "C_tier_tests.tsv"), sep="\t")
        c3 = C[(C.scope == "tier3_v_tier0") &
               (C.metric == "phylop_mean")].iloc[0]
        rows = [
            "| Evo-ML ablation (28A) | +phyloP {} macro-F1 ({:+.4f} vs F0 "
            "0.336; phyloP-replaces-phast {:.4f}); F0 gate passed | "
            "results/evo_ext/A_ablation.tsv |".format(
                "improved" if f2.delta_vs_F0 > 0 else "did not improve",
                f2.delta_vs_F0, f3.macro_f1),
            "| Broad vs specific (28B) | GTEx broad v specific lnc: phylop "
            "p={:.2e} (median {} v {}) | results/evo_ext/B_*.tsv |".format(
                bl.pval, bl.median_broad, bl.median_specific),
            "| Cancer DE tiers (28C) | tier3(up>=5) v tier0 phylop p={:.2e}; "
            "18-cancer per-gene DE, 13d spot-check passed | "
            "results/evo_ext/C_*.tsv |".format(c3.pval),
        ]
        path = os.path.join(WORKDIR, "REPORT.md")
        with open(path) as fh:
            text = fh.read()
        with open(safe_write_path(path), "w") as fh:
            fh.write(text.rstrip("\n") + "\n" + "\n".join(rows) + "\n")
        st.note("REPORT.md: appended {} evolution-extension rows".format(
            len(rows)))


if __name__ == "__main__":
    mode = (sys.argv[1] if len(sys.argv) > 1 else "").upper()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    if mode == "A":
        mode_a()
    elif mode == "B":
        mode_b()
    elif mode == "C":
        mode_c()
    elif mode == "R":
        mode_r()
    else:
        raise SystemExit("usage: 28_evo_exts.py A|B|C|R")
