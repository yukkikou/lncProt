#!/usr/bin/env python
"""13d_tcga_disease.py — disease dimension (TCGA), PARALLEL per cancer type.

Modes:
  one <CANCER>   -> compute per-cancer stats, write shard
                    data/tissue_disease/shards/13d_<CANCER>.tsv
  merge          -> concatenate shards, global class-level summary,
                    LTR-derepression ranking

Per cancer (tumor n>=20; normal any n):
  - per-CLASS differential expression tumor vs normal (Mann-Whitney per gene,
    BH within cancer; aggregated: frac up/down, median delta log2TPM, class-level
    MWU of per-gene deltas vs PCG background)
  - LTR-derived lncRNA derepression (median delta of LTR-lnc vs other lnc,
    enrichment of top-up genes among LTR-lnc; Fisher)
  - sample accounting row (n tumor, n normal, n genes tested)
Outputs -> results/tissue_disease/tables/tcga_class_de_per_cancer.tsv,
           tcga_ltr_derepression.tsv, tcga_class_stage.tsv (stage if n>=10/stage)
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

TD = os.path.join(WORKDIR, "data", "tissue_disease")
SHARDS = os.path.join(TD, "shards")
OUT = os.path.join(WORKDIR, "results", "tissue_disease", "tables")
MIN_TUMOR = 20
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


def load():
    d = np.load(os.path.join(TD, "toil_class_expr.npz"))["M"]
    genes = pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t", index_col=0)
    samp = pd.read_csv(os.path.join(TD, "toil_samples.tsv"), sep="\t")
    return d, genes, samp


def one(cancer):
    with step(f"13d_{cancer}") as st:
        d, genes, samp = load()
        tum = samp[(samp.cancer == cancer) & (samp.condition == "tumor")].index
        nor = samp[(samp.cancer == cancer) & (samp.condition == "normal")].index
        if len(tum) < MIN_TUMOR:
            st.note(f"skip: only {len(tum)} tumors")
            return
        d = np.nan_to_num(d)
        T = d[:, tum.values].mean(axis=1)
        N = d[:, nor.values].mean(axis=1) if len(nor) >= 5 else None
        detect = T > 0 if N is None else (T > 0) | (N > 0)
        delta = (T - N) if N is not None else np.full(len(genes), np.nan)
        delta[~detect] = np.nan

        # per-gene MWU + BH within this cancer (only when normal comparator)
        pvals = np.full(len(genes), np.nan)
        if N is not None:
            for i in range(len(genes)):
                a = d[i, tum.values]
                b = d[i, nor.values]
                if a.std() == 0 and b.std() == 0:
                    continue
                _, p = mannwhitneyu(a, b, alternative="two-sided")
                pvals[i] = p
        fdr = bh_fdr(pvals[np.isfinite(pvals)])
        fdr_full = np.full(len(genes), np.nan)
        fdr_full[np.isfinite(pvals)] = fdr

        rows, ltr_rows = [], []
        lnc_mask = genes.cls.isin(CLASSES).values
        for c in CLASSES + ["PCG"]:
            idx = (genes.cls == "PCG").values if c == "PCG" else \
                  (genes.cls == c).values
            if N is not None:
                up = (delta > 0.5) & (fdr_full < 0.05)
                dn = (delta < -0.5) & (fdr_full < 0.05)
                den = (idx & detect).sum()
                fu, fd = float((up & idx).sum() / max(den, 1)), \
                         float((dn & idx).sum() / max(den, 1))
                med = float(np.nanmedian(delta[idx]))
            else:
                fu = fd = med = np.nan
            rows.append(dict(cancer=cancer, cls=c, n_genes=int(idx.sum()),
                             n_tumor=len(tum), n_normal=len(nor),
                             median_delta=med, frac_up=fu, frac_down=fd))
            # LTR derepression stats once (lnc only)
            if c == "intergenic":
                ltr = lnc_mask & (genes.tss_in_ltr == 1).values
                other = lnc_mask & ~(genes.tss_in_ltr == 1).values
                if N is not None and ltr.sum() >= 10:
                    m_l = float(np.nanmedian(delta[ltr]))
                    m_o = float(np.nanmedian(delta[other]))
                    up = (delta > 0.5) & (fdr_full < 0.05)
                    a11 = int((up & ltr).sum())
                    a12 = int((~up & ltr).sum())
                    a21 = int((up & other).sum())
                    a22 = int((~up & other).sum())
                    orr, p = fisher_exact([[a11, a12], [a21, a22]])
                    ltr_rows.append(dict(cancer=cancer, n_ltr_lnc=int(ltr.sum()),
                                         median_delta_LTR=m_l, median_delta_other=m_o,
                                         diff=m_l - m_o, fisher_OR=orr, fisher_p=p))
        os.makedirs(SHARDS, exist_ok=True)
        pd.DataFrame(rows).to_csv(safe_write_path(
            os.path.join(SHARDS, f"13d_{cancer}.tsv")), sep="\t", index=False)
        if ltr_rows:
            pd.DataFrame(ltr_rows).to_csv(safe_write_path(
                os.path.join(SHARDS, f"13d_ltr_{cancer}.tsv")), sep="\t",
                index=False)
        st.note(f"tumor={len(tum)} normal={len(nor)}; classes={len(rows)}")


def merge():
    with step("13d_merge") as st:
        os.makedirs(OUT, exist_ok=True)
        shards = sorted(f for f in os.listdir(SHARDS) if f.startswith("13d_")
                        and not f.startswith("13d_ltr_") and f.endswith(".tsv"))
        de = pd.concat([pd.read_csv(os.path.join(SHARDS, f), sep="\t")
                        for f in shards], ignore_index=True)
        de.to_csv(safe_write_path(os.path.join(
            OUT, "tcga_class_de_per_cancer.tsv")), sep="\t", index=False)
        st.note(f"cancers={de.cancer.nunique()} rows={len(de)}")

        lshards = sorted(f for f in os.listdir(SHARDS) if f.startswith("13d_ltr_"))
        if lshards:
            lt = pd.concat([pd.read_csv(os.path.join(SHARDS, f), sep="\t")
                            for f in lshards], ignore_index=True)
            lt = lt.sort_values("diff", ascending=False)
            lt.to_csv(safe_write_path(os.path.join(
                OUT, "tcga_ltr_derepression.tsv")), sep="\t", index=False)
            st.note("top LTR derepression: " + "; ".join(
                f"{r.cancer}: d={r.diff:.2f}" for r in lt.head(5).itertuples()))
        # class summary across cancers
        sm = (de[de.cls != "PCG"].groupby("cls")
              .agg(n_cancers=("cancer", "nunique"),
                   median_delta=("median_delta", "median"),
                   frac_up=("frac_up", "median")).reset_index())
        st.summary_md("results/tissue_disease/step_summaries/13d_tcga_disease.md",
                      "13d TCGA disease dimension",
                      "## Class DE across cancers\n\n" + df_to_md(sm, ".3f")
                      + ("\n\n## LTR-lncRNA derepression (top)\n\n"
                         + df_to_md(lt.head(15), ".3f") if lshards else ""))


def stage():
    """Stage analysis via LOCAL clinical_simple.tsv (case_submitter_id +
    AJCC pathologic stage). NOTE: stage_anno.txt file UUIDs belong to a
    different file batch than gdc_sample_sheet (zero overlap, verified), so
    stage is joined at the patient level from clinical tables instead."""
    with step("13d_stage") as st:
        clin = pd.read_csv(os.path.join(TD, "clinical_simple.tsv"), sep="\t",
                           low_memory=False)
        clin["patient12"] = clin["case_submitter_id"].astype(str).str[:12]
        clin["stage"] = clin["ajcc_pathologic_stage"].fillna(
            clin["ajcc_clinical_stage"])
        clin = clin.dropna(subset=["stage"])
        smap = clin.drop_duplicates("patient12").set_index("patient12")["stage"]
        d, genes, samp = load()
        d = np.nan_to_num(d)
        rows = []
        for cancer, sub in samp[(samp.condition == "tumor") &
                                samp.cancer.notna()].groupby("cancer"):
            stg = sub["sample"].str[:12].map(smap)
            sub = sub.assign(stage=stg.values)
            early = sub[sub.stage.str.contains(
                r"Stage I\b|Stage II\b|Stage 0", na=False, regex=True)].index
            adv = sub[sub.stage.str.contains(
                "Stage III|Stage IV", na=False)].index
            if len(early) < 10 or len(adv) < 10:
                continue
            dE = d[:, sub.loc[early].index.values].mean(axis=1)
            dA = d[:, sub.loc[adv].index.values].mean(axis=1)
            detect = (dE > 0) | (dA > 0)
            dd = np.where(detect, dA - dE, np.nan)
            for c in CLASSES:
                idx = (genes.cls == c).values
                rows.append(dict(cancer=cancer, cls=c, n_early=len(early),
                                 n_advanced=len(adv),
                                 median_delta_advanced=float(np.nanmedian(dd[idx]))))
        if rows:
            t = pd.DataFrame(rows)
            t.to_csv(safe_write_path(os.path.join(OUT, "tcga_class_stage.tsv")),
                     sep="\t", index=False)
            st.note(f"stage rows={len(t)} across {t.cancer.nunique()} cancers")
        else:
            st.note("no cancer met stage thresholds (n>=10 early & advanced)")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "merge"
    if mode == "one":
        one(sys.argv[2])
    elif mode == "stage":
        stage()
    else:
        merge()
