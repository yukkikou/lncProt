#!/usr/bin/env python
"""33_cancer_bridge.py — bridge to Sun et al. Sci Adv 2026 definition.

Reproduce a pan-cancer-associated lncRNA (PC-lncRNA) definition on top of
our birth classes:
  - per-cancer per-gene DE (tumour vs matched normal; 13d logic), saved
  - per-cancer per-gene survival association (vectorised log-rank)
  - PC_like = DE in >=half of cancers AND survival P<0.05 in >=1 cancer
  - class composition + logistic regression testing whether birth class
    adds beyond expression breadth (CAGE breadth, tau, max TPM)
Outputs -> results/cancer_bridge/*.tsv + step summary.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist
from scipy.stats import fisher_exact, mannwhitneyu
import statsmodels.api as sm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
TD = os.path.join(WORKDIR, "data", "tissue_disease")
OUT = os.path.join(RES, "cancer_bridge")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
MIN_TUMOR, MIN_NORMAL = 20, 5


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def load_survival():
    c = pd.read_csv(os.path.join(TD, "clinical_simple.tsv"), sep="\t",
                    low_memory=False)
    c["patient"] = c["case_submitter_id"].astype(str).str[:12]
    dead = c["vital_status"].astype(str).str.lower().eq("dead")
    t = pd.to_numeric(c["days_to_death"], errors="coerce").fillna(
        pd.to_numeric(c["days_to_last_follow_up"], errors="coerce"))
    out = pd.DataFrame({"patient": c["patient"], "time": t,
                        "event": dead.astype(int)})
    out = out.dropna(subset=["time"])
    return out[out.time > 0].drop_duplicates("patient").set_index("patient")


def logrank_matrix(G, time, event):
    """Vectorised two-sample log-rank across rows of G (1/0 group)."""
    order = np.argsort(time)
    t = np.asarray(time)[order]
    e = np.asarray(event)[order]
    Gs = G[:, order]
    total1 = Gs.sum(axis=1)
    cs = np.cumsum(Gs, axis=1)
    uniq, first = np.unique(t, return_index=True)
    O1 = np.zeros(Gs.shape[0], dtype=float)
    E1 = np.zeros(Gs.shape[0], dtype=float)
    V = np.zeros(Gs.shape[0], dtype=float)
    n_total = len(t)
    for tk, i0 in zip(uniq, first):
        idx = np.where(t == tk)[0]
        d = int(e[idx].sum())
        if d == 0:
            continue
        n = n_total - i0
        n1 = total1 if i0 == 0 else total1 - cs[:, i0 - 1]
        d1 = Gs[:, idx].sum(axis=1)
        frac = n1 / n
        O1 += d1
        E1 += d * frac
        if n > 1:
            V += d * frac * (1 - frac) * (n - d) / (n - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        stat = (O1 - E1) ** 2 / V
        p = chi2_dist.sf(stat, 1)
    p[~np.isfinite(p)] = np.nan
    return p, O1, E1


def main():
    with step("33_cancer_bridge") as st:
        d = np.nan_to_num(np.load(os.path.join(TD, "toil_class_expr.npz"))
                          ["M"])
        genes = pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t",
                            index_col=0)
        samp = pd.read_csv(os.path.join(TD, "toil_samples.tsv"), sep="\t")
        lnc_mask = genes.cls.isin(CLASSES).values
        lnc_idx = np.where(lnc_mask)[0]
        lnc_ids = genes.index.values[lnc_idx]
        lnc_vids = genes["gene_id"].values[lnc_idx]   # versioned ids
        n_lnc = len(lnc_idx)
        D = d[lnc_idx]                      # one copy, reused below

        cancers = [c for c, sub in samp[samp.cancer.notna()].groupby(
            "cancer") if ((sub.condition == "tumor").sum() >= MIN_TUMOR
                          and (sub.condition == "normal").sum() >= MIN_NORMAL)]
        st.note("cancers with matched normals: {}".format(len(cancers)))

        delta_m = np.full((n_lnc, len(cancers)), np.nan)
        fdr_m = np.full((n_lnc, len(cancers)), np.nan)
        pv_m = np.full((n_lnc, len(cancers)), np.nan)
        for j, cancer in enumerate(cancers):
            tum = samp[(samp.cancer == cancer) &
                       (samp.condition == "tumor")].index.values
            nor = samp[(samp.cancer == cancer) &
                       (samp.condition == "normal")].index.values
            T = D[:, tum].mean(axis=1)
            N = D[:, nor].mean(axis=1)
            detect = (T > 0) | (N > 0)
            delta = np.where(detect, T - N, np.nan)
            pv = np.full(n_lnc, np.nan)
            for i in range(n_lnc):
                a = D[i, tum]
                b = D[i, nor]
                if a.std() == 0 and b.std() == 0:
                    continue
                from scipy.stats import mannwhitneyu
                _, p = mannwhitneyu(a, b, alternative="two-sided")
                pv[i] = p
            fdr = np.full(n_lnc, np.nan)
            ok = np.isfinite(pv)
            fdr[ok] = bh_fdr(pv[ok])
            delta_m[:, j] = delta
            fdr_m[:, j] = fdr
            pv_m[:, j] = pv
            st.note("DE {}: tumour={} normal={}".format(
                cancer, len(tum), len(nor)))

        de_up = (delta_m > 0.5) & (fdr_m < 0.05)
        de_dn = (delta_m < -0.5) & (fdr_m < 0.05)
        n_up = np.nansum(de_up, axis=1)
        n_dn = np.nansum(de_dn, axis=1)
        de_any = de_up | de_dn
        n_de = np.nansum(de_any, axis=1)

        # ---- survival per gene per cancer (tumours with survival) ---------
        surv = load_survival()
        samp["patient"] = samp["sample"].astype(str).str[:12]
        samp = samp.join(surv, on="patient")
        surv_p = np.full((n_lnc, len(cancers)), np.nan)
        for j, cancer in enumerate(cancers):
            sub = samp[(samp.cancer == cancer) &
                       (samp.condition == "tumor") & samp.time.notna()]
            if len(sub) < 40 or int(sub.event.sum()) < 15:
                continue
            E = D[:, sub.index.values]
            G = (E > np.median(E, axis=1, keepdims=True)).astype(float)
            p, _, _ = logrank_matrix(G, sub.time.values, sub.event.values)
            surv_p[:, j] = p
            st.note("surv {}: n={} events={} p<0.05 genes={}".format(
                cancer, len(sub), int(sub.event.sum()),
                int(np.nansum(p < 0.05))))
        n_surv_sig = np.nansum(surv_p < 0.05, axis=1)

        half = int(np.ceil(len(cancers) / 2))
        pc_like = (n_de >= half) & (n_surv_sig >= 1)
        pc_like_strict = (n_de >= half) & (n_surv_sig >= 2)

        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        act = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                          index_col=0, low_memory=False)[
            ["breadth", "tau", "max_tpm", "mean_tpm", "active"]]
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)[
            ["phylop_mean", "phast_mean", "gc", "tss_in_ltr"]]
        per = pd.DataFrame({
            "gene_id": lnc_vids,
            "base_id": lnc_ids,
            "cls": genes.loc[lnc_ids, "cls"].values,
            "n_de": n_de, "n_up": n_up, "n_dn": n_dn,
            "n_surv_sig": n_surv_sig,
            "pc_like": pc_like.astype(int),
            "pc_like_strict": pc_like_strict.astype(int)},
            index=lnc_vids)
        per = per.join(act).join(ph)
        per.to_csv(safe_write_path(os.path.join(
            OUT, "C_bridge_pergene.tsv")), sep="\t")

        rows = []
        tab = pd.crosstab(per.cls, per.pc_like).reindex(
            columns=[0, 1], fill_value=0)
        from scipy.stats import chi2_contingency
        try:
            chi2, pchi, _, _ = chi2_contingency(tab.values)
        except Exception:
            chi2, pchi = np.nan, np.nan
        rows.append(dict(analysis="composition", contrast="PC_like x class",
                         test="chi2", stat=chi2, pval=pchi))
        for c in CLASSES:
            a = int(tab.loc[c, 1]) if c in tab.index else 0
            b = int(tab.loc[c, 0]) if c in tab.index else 0
            rest_a = int(per[(per.cls != c) & (per.pc_like == 1)].shape[0])
            rest_b = int(per[(per.cls != c) & (per.pc_like == 0)].shape[0])
            orr, pf = fisher_exact([[a, b], [rest_a, rest_b]])
            rows.append(dict(analysis="composition", contrast="OR_" + c,
                             test="fisher", stat=orr, pval=pf))
        # direction among PC-like
        pcg_like = per[per.pc_like == 1]
        rows.append(dict(analysis="direction", contrast="PC_like_up_vs_down",
                         test="count", stat=float((pcg_like.n_up >
                                                   pcg_like.n_dn).mean()),
                         pval=np.nan))
        # logistic regression: class + breadth
        X = per.copy()
        X["log_breadth"] = np.log1p(X.breadth)
        X["log_tpm"] = np.log1p(X.max_tpm)
        X = X.dropna(subset=["log_breadth", "tau", "log_tpm", "pc_like"])
        dummies = pd.get_dummies(X.cls, prefix="cls", drop_first=True)
        Xb = pd.concat([X[["log_breadth", "tau", "log_tpm"]], dummies],
                       axis=1).astype(float)
        y = X.pc_like.astype(float)
        try:
            m_full = sm.Logit(y, sm.add_constant(Xb)).fit(disp=0)
            m_red = sm.Logit(y, sm.add_constant(
                X[["log_breadth", "tau", "log_tpm"]].astype(float))
            ).fit(disp=0)
            lr = 2 * (m_full.llf - m_red.llf)
            df = Xb.shape[1] - 3
            p_lr = chi2_dist.sf(lr, df)
            rows.append(dict(analysis="logistic", contrast="class_beyond_breadth",
                             test="LR", stat=float(lr), pval=float(p_lr)))
            for name in m_full.params.index:
                if name == "const":
                    continue
                rows.append(dict(analysis="logistic_coef", contrast=name,
                                 test="coef", stat=float(m_full.params[name]),
                                 pval=float(m_full.pvalues[name])))
        except Exception as e:
            st.note("logistic failed: {}".format(e))
        T = pd.DataFrame(rows)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "C_bridge_tests.tsv")), sep="\t", index=False)
        st.note("PC_like n={} ({}% of lnc); up>down frac={:.2f}; class LR "
                "p={}".format(int(per.pc_like.sum()),
                              100 * per.pc_like.mean(),
                              float((pcg_like.n_up > pcg_like.n_dn).mean()),
                              T[(T.analysis == "logistic")].pval.iloc[0]
                              if (T.analysis == "logistic").any() else np.nan))
        st.summary_md("results/cancer_bridge/step_summaries/33_bridge.md",
                      "33 PC-lncRNA definition x birth class",
                      "## Class composition of PC-like lncRNAs\n\n"
                      + df_to_md(pd.crosstab(per.cls,
                                             per.pc_like).reset_index(), "d")
                      + "\n\n## Tests\n\n" + df_to_md(T, ".4g", max_rows=40)
                      + "\n\nPC_like = DE in >= {} / {} cancers AND "
                        "log-rank P<0.05 in >=1 cancer.".format(
                            half, len(cancers)))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    main()
