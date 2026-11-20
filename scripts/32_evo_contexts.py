#!/usr/bin/env python
"""32_evo_contexts.py — evolutionary/clinical context analyses inspired by
Sun et al. Sci Adv 2026 (rapid evolution of cancer-associated lncRNAs).

Modes:
  A  TCGA survival by birth class (per-cancer median split log-rank +
     continuous Cox per SD; Stouffer meta across cancers).
  B  Recent selection at promoters: 1000G Tajima's D (CEU/CHB/YRI) mapped
     to -500/+200 promoter windows; class contrasts + LTR subset.
  C  Transcript isoform number / length by birth class (GENCODE).
  D  Partner-PCG age proxy (conservation) for divergent loci.
Outputs -> results/evo_contexts/*.tsv + step summaries.
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.duration.survfunc import survdiff

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
TD = os.path.join(WORKDIR, "data", "tissue_disease")
OUT = os.path.join(RES, "evo_contexts")
TAD = "/media/share/resource/database/1000_Genomes_Selection"
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


def load_survival():
    c = pd.read_csv(os.path.join(TD, "clinical_simple.tsv"), sep="\t",
                    low_memory=False)
    c["patient"] = c["case_submitter_id"].astype(str).str[:12]
    dead = c["vital_status"].astype(str).str.lower().eq("dead")
    t = pd.to_numeric(c["days_to_death"], errors="coerce")
    t = t.fillna(pd.to_numeric(c["days_to_last_follow_up"],
                               errors="coerce"))
    out = pd.DataFrame({"patient": c["patient"], "time": t,
                        "event": dead.astype(int)})
    out = out.dropna(subset=["time"])
    out = out[out.time > 0].drop_duplicates("patient")
    return out.set_index("patient")


# ------------------------------------------------------------------ mode A --
def mode_a():
    with step("32_ctx_A") as st:
        d = np.nan_to_num(np.load(os.path.join(TD, "toil_class_expr.npz"))
                          ["M"])
        genes = pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t",
                            index_col=0)
        samp = pd.read_csv(os.path.join(TD, "toil_samples.tsv"), sep="\t")
        surv = load_survival()
        samp["patient"] = samp["sample"].astype(str).str[:12]
        samp = samp.join(surv, on="patient")
        tum = samp[(samp.cohort == "TCGA") & (samp.condition == "tumor")
                   & samp.time.notna()]
        lnc = genes.cls.isin(CLASSES).values
        ltr = lnc & (genes.tss_in_ltr == 1).fillna(False).values
        oth = lnc & ~ltr
        scores = {}
        for c in CLASSES:
            scores[c] = d[(genes.cls == c).values].mean(axis=0)
        scores["LTR_lnc"] = d[ltr].mean(axis=0)
        scores["other_lnc"] = d[oth].mean(axis=0)
        sc = pd.DataFrame(scores, index=samp.index)
        rows = []
        for cancer, sub in tum.groupby("cancer"):
            if len(sub) < 40 or int(sub.event.sum()) < 15:
                continue
            for feat in scores:
                v = sc.loc[sub.index, feat]
                if v.std() == 0:
                    continue
                # per-SD Cox via statsmodels PHReg
                try:
                    from statsmodels.duration.hazard_regression import PHReg
                    z = (v - v.mean()) / v.std()
                    m = PHReg(sub.time.values, z.values,
                              status=sub.event.values).fit(disp=0)
                    hr, pcox = float(np.exp(m.params[0])), float(m.pvalues[0])
                except Exception:
                    hr, pcox = np.nan, np.nan
                # median-split log-rank
                g = (v > v.median()).astype(int).values
                try:
                    stat, plr = survdiff(sub.time.values,
                                         sub.event.values, g)
                    plr = float(plr)
                except Exception:
                    plr = np.nan
                rows.append(dict(cancer=cancer, feature=feat, n=len(sub),
                                 events=int(sub.event.sum()), HR_perSD=hr,
                                 p_cox=pcox, p_logrank=plr,
                                 direction="risk" if hr > 1 else "protective"))
        T = pd.DataFrame(rows)
        T["fdr_logrank"] = np.nan
        T["fdr_cox"] = np.nan
        T["fdr_logrank"] = bh_fdr(T.p_logrank.fillna(1).values)
        T["fdr_cox"] = bh_fdr(T.p_cox.fillna(1).values)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "A_survival_per_cancer.tsv")), sep="\t", index=False)
        # Stouffer meta per feature (sign from Cox HR)
        meta = []
        for feat, sub in T.groupby("feature"):
            p = sub.p_cox.dropna().values
            hr = sub.HR_perSD.dropna().values
            if len(p) == 0:
                continue
            z = np.sign(np.log(hr)) * -np.sqrt(2) * np.log(p / 2)
            zmeta = z.sum() / np.sqrt(len(z))
            from scipy.stats import norm
            pmeta = 2 * norm.sf(abs(zmeta))
            meta.append(dict(feature=feat, n_cancers=len(z),
                             n_sig_fdr=float((sub.fdr_cox < 0.05).sum()),
                             n_risk=float((sub.HR_perSD > 1).sum()),
                             n_protective=float((sub.HR_perSD < 1).sum()),
                             z_meta=float(zmeta), p_meta=float(pmeta)))
        M = pd.DataFrame(meta)
        M.to_csv(safe_write_path(os.path.join(
            OUT, "A_survival_meta.tsv")), sep="\t", index=False)
        st.note("A survival: cancers={}; meta: ".format(
            T.cancer.nunique()) + "; ".join(
            "{} z={:.2f} p={:.1e}".format(r.feature, r.z_meta, r.p_meta)
            for r in M.itertuples()))
        st.summary_md("results/evo_contexts/step_summaries/32_a.md",
                      "32A TCGA survival by birth class",
                      "## Meta across cancers (Cox per SD)\n\n"
                      + df_to_md(M, ".3f")
                      + "\n\n## Per-cancer top rows\n\n"
                      + df_to_md(T.sort_values("p_cox").head(25), ".3f"))


# ------------------------------------------------------------------ mode B --
def load_tajima(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        t = pd.read_csv(fh, sep=r"\s+", header=0)
    t.columns = [c.strip() for c in t.columns]
    t["chrom"] = t["chromosome"].astype(str).str.replace("chr", "",
                                                         regex=False)
    return t[["chrom", "start", "end", "score", "pvalue"]]


def promoter_tajima(prom, tad):
    """Overlap-weighted mean Tajima's D per promoter."""
    out = np.full(len(prom), np.nan)
    out_p = np.full(len(prom), np.nan)
    for chrom, idx in prom.groupby("chrom").groups.items():
        sub = tad[tad.chrom == chrom.replace("chr", "")]
        if len(sub) == 0:
            continue
        starts = sub.start.values
        ends = sub.end.values
        score = sub.score.values
        pv = sub.pvalue.values
        p = prom.loc[idx]
        for j, (s, e) in enumerate(zip(p.start.values, p.end.values)):
            i0 = np.searchsorted(ends, s, side="right")
            i1 = np.searchsorted(starts, e, side="left")
            if i1 <= i0:
                continue
            ov = (np.minimum(ends[i0:i1], e) -
                  np.maximum(starts[i0:i1], s)).clip(min=0)
            w = ov.sum()
            if w <= 0:
                continue
            out[idx[j]] = float((score[i0:i1] * ov).sum() / w)
            out_p[idx[j]] = float((pv[i0:i1] * ov).sum() / w)
    return out, out_p


def mode_b():
    with step("32_ctx_B") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        prom = pd.read_csv(os.path.join(RES, "promoters_feat.bed"), sep="\t",
                           names=["chrom", "start", "end", "gene_id", "x",
                                  "strand"])
        cols = {}
        for pop in ["CEU", "CHB", "YRI"]:
            f = os.path.join(TAD, "TajimasD_{}.whole_genome.pvalues".format(
                pop))
            if not os.path.exists(f):
                f += ".gz"
            t = load_tajima(f)
            m, mp = promoter_tajima(prom, t)
            cols[pop + "_tajD"] = m
            cols[pop + "_p"] = mp
            st.note("{} windows={}".format(pop, len(t)))
        D = pd.DataFrame(cols, index=prom.gene_id.values)
        D["tajD_mean"] = D[[c + "_tajD" for c in ["CEU", "CHB", "YRI"]]] \
            .mean(axis=1)
        D["p_mean"] = D[[c + "_p" for c in ["CEU", "CHB", "YRI"]]].mean(axis=1)
        m = meta.join(D, how="left")
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"),
                         sep="\t", index_col=0, low_memory=False)[
            ["tss_in_ltr"]]
        m = m.join(ph[["phylop_mean"]]).join(te)
        m[["tajD_mean", "p_mean"]].to_csv(safe_write_path(os.path.join(
            OUT, "B_tajima_per_gene.tsv")), sep="\t")
        rows = []
        pcg = m[m.group == "pcg"]
        for feat in ["tajD_mean", "p_mean"]:
            x = pcg[feat].dropna().values
            for cls in CLASSES + ["pseudogene", "transcribed_pseudogene"]:
                y = m.loc[m.cls == cls, feat].dropna().values
                if len(y) < 20:
                    continue
                dlt, p = mwu(y, x)
                rows.append(dict(feature=feat, cls=cls, n=len(y),
                                 median_class=float(np.median(y)),
                                 median_PCG=float(np.median(x)),
                                 cliffs_delta=dlt, pval=p))
        # LTR vs non-LTR within intergenic
        ig = m[m.cls == "intergenic"]
        a = ig.loc[ig.tss_in_ltr == 1, "tajD_mean"].dropna().values
        b = ig.loc[ig.tss_in_ltr != 1, "tajD_mean"].dropna().values
        if len(a) >= 20 and len(b) >= 20:
            dlt, p = mwu(a, b)
            rows.append(dict(feature="tajD_mean",
                             cls="intergenic_LTR_vs_nonLTR", n=len(a) + len(b),
                             median_class=float(np.median(a)),
                             median_PCG=float(np.median(b)),
                             cliffs_delta=dlt, pval=p))
        T = pd.DataFrame(rows)
        T["fdr"] = bh_fdr(T.pval.values)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "B_tajima_tests.tsv")), sep="\t", index=False)
        ok = m[["tajD_mean", "phylop_mean"]].dropna()
        r, pr = spearmanr(ok.tajD_mean, ok.phylop_mean)
        st.note("B Tajima: TajD mean={:.3f}; corr(phyloP)=".format(
            m.tajD_mean.mean()) + "rho={:.3f} p={:.1e}; LTR vs nonLTR "
            "d={:.3f} p={:.1e}".format(
                r, pr, T.iloc[-1].cliffs_delta, T.iloc[-1].pval))
        st.summary_md("results/evo_contexts/step_summaries/32_b.md",
                      "32B 1000G Tajima's D at promoters",
                      "## Class vs PCG\n\n" + df_to_md(T, ".3f")
                      + "\n\nCorrelation with phyloP: rho={:.3f}, "
                        "P={:.2e}".format(r, pr))


# ------------------------------------------------------------------ mode C --
def mode_c():
    with step("32_ctx_C") as st:
        tx = pd.read_csv(os.path.join(WORKDIR, "data", "genes",
                                      "transcripts.tsv"), sep="\t")
        agg = tx.groupby("gene_id").agg(
            n_tx=("tx_id", "nunique"),
            max_len=("length", "max"),
            med_len=("length", "median"),
            max_exons=("n_exons", "max"),
            has_mane=("mane", lambda s: int(pd.to_numeric(
                s, errors="coerce").fillna(0).max() > 0)))
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"),
                         sep="\t", index_col=0, low_memory=False)[
            ["tss_in_ltr"]]
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)[
            ["phylop_mean"]]
        m = meta[["cls", "group"]].join(agg).join(te).join(ph)
        rows = []
        pcg = m[m.group == "pcg"]
        for feat in ["n_tx", "max_len", "med_len", "max_exons"]:
            x = pcg[feat].dropna().values
            for cls in CLASSES + ["pseudogene", "transcribed_pseudogene"]:
                y = m.loc[m.cls == cls, feat].dropna().values
                if len(y) < 20:
                    continue
                dlt, p = mwu(y, x)
                rows.append(dict(feature=feat, cls=cls, n=len(y),
                                 median_class=float(np.median(y)),
                                 median_PCG=float(np.median(x)),
                                 cliffs_delta=dlt, pval=p))
        ig = m[m.cls == "intergenic"]
        for feat in ["n_tx", "max_len", "max_exons"]:
            a = ig.loc[ig.tss_in_ltr == 1, feat].dropna().values
            b = ig.loc[ig.tss_in_ltr != 1, feat].dropna().values
            if len(a) >= 20 and len(b) >= 20:
                dlt, p = mwu(a, b)
                rows.append(dict(feature=feat,
                                 cls="intergenic_LTR_vs_nonLTR",
                                 n=len(a) + len(b),
                                 median_class=float(np.median(a)),
                                 median_PCG=float(np.median(b)),
                                 cliffs_delta=dlt, pval=p))
        T = pd.DataFrame(rows)
        T["fdr"] = bh_fdr(T.pval.values)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "C_isoform_tests.tsv")), sep="\t", index=False)
        st.note("C isoforms: lnc n_tx median={}; PCG={}; LTR-vs-nonLTR "
                "n_tx d={:.3f}".format(
                    m[m.group == "lnc"].n_tx.median(),
                    pcg.n_tx.median(),
                    T[(T.cls == "intergenic_LTR_vs_nonLTR") &
                      (T.feature == "n_tx")].cliffs_delta.iloc[0]))
        st.summary_md("results/evo_contexts/step_summaries/32_c.md",
                      "32C isoform number / transcript length",
                      "## Class vs PCG\n\n" + df_to_md(T, ".3f"))


# ------------------------------------------------------------------ mode D --
def mode_d():
    with step("32_ctx_D") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)
        m = meta.join(ph[["phylop_mean", "phast_mean"]])
        partners = set(meta.loc[meta.partner_pcg.notna() &
                                meta.partner_pcg.astype(str).ne(""),
                                "partner_pcg"].astype(str))
        pcg = m[m.group == "pcg"].copy()
        pcg["is_partner"] = pcg.index.isin(partners)
        rows = []
        for feat in ["phylop_mean", "phast_mean"]:
            a = pcg.loc[pcg.is_partner, feat].dropna().values
            b = pcg.loc[~pcg.is_partner, feat].dropna().values
            dlt, p = mwu(a, b)
            rows.append(dict(feature=feat, contrast="partner_vs_nonpartner_PCG",
                             n_partner=len(a), n_other=len(b),
                             median_partner=float(np.median(a)),
                             median_other=float(np.median(b)),
                             cliffs_delta=dlt, pval=p))
        # within divergent: own promoter vs partner PCG conservation
        div = m[(m.cls == "divergent") & m.partner_pcg.notna()].copy()
        div["partner"] = div.partner_pcg.astype(str)
        pv = pcg[["phylop_mean", "phast_mean"]].rename(
            columns={"phylop_mean": "partner_phylop",
                     "phast_mean": "partner_phast"})
        div = div.join(pv, on="partner")
        for a, b in [("phylop_mean", "partner_phylop"),
                     ("phast_mean", "partner_phast")]:
            ok = div[[a, b]].dropna()
            if len(ok) >= 30:
                r, pr = spearmanr(ok[a], ok[b])
                rows.append(dict(feature=a, contrast="div_vs_partner_rho",
                                 n_partner=len(ok), n_other=0,
                                 median_partner=np.nan, median_other=np.nan,
                                 spearman_rho=float(r), pval=float(pr)))
        T = pd.DataFrame(rows)
        T["fdr"] = bh_fdr(T.pval.values)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "D_partner_age_proxy.tsv")), sep="\t", index=False)
        st.note("D partner age proxy: " + "; ".join(
            "{} {} d={:.3f} p={:.1e}".format(r.contrast, r.feature,
                                             getattr(r, "cliffs_delta", np.nan),
                                             r.pval)
            for r in T.itertuples()))
        st.summary_md("results/evo_contexts/step_summaries/32_d.md",
                      "32D partner-PCG conservation (age proxy)",
                      "## Partner vs non-partner PCGs\n\n"
                      + df_to_md(T, ".3f"))


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
    elif mode == "D":
        mode_d()
    else:
        raise SystemExit("usage: 32_evo_contexts.py A|B|C|D")
