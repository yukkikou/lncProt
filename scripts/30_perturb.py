#!/usr/bin/env python
"""30_perturb.py — perturbation validation of two Discussion predictions.

Prediction 1 (DUX4): experimental DUX4 activation should selectively induce
  LTR-derived lncRNA transcription.
  - GSE85461: doxycycline-inducible DUX4 in myoblasts (provided edgeR logFC).
  - GSE56787: FSHD vs control cultured myoblasts (endogenous DUX4; counts).
Prediction 2 (hypomethylation): LTR-lncRNA activation should track retroviral
  chromatin opening under DNA hypomethylation.
  - GSE337536: MB2141 melanoma untreated vs decitabine (2x2, counts).
  - GSE305426: MPNST lines DMSO vs decitabine/azacitidine (featureCounts).

Per dataset: per-gene response (log2FC), LTR-derived lnc vs other lnc (MWU +
Fisher on up fraction), PCG comparator row. Outputs -> results/perturb/.
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
TD = os.path.join(WORKDIR, "data", "tissue_disease")
PD = os.path.join(WORKDIR, "data", "perturb")
OUT = os.path.join(RES, "perturb")


def load_classes():
    tg = pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t",
                     index_col=0)
    lab = tg[["gene_id", "cls", "tss_in_ltr"]].copy()
    lab["group3"] = np.where(
        ~lab.cls.isin(["divergent", "antisense", "intronic",
                       "sense_overlapping", "intergenic"]), "other",
        np.where(lab.tss_in_ltr.fillna(0) == 1, "LTR_lnc", "lnc_other"))
    lab.loc[lab.cls == "PCG", "group3"] = "PCG"
    return lab.set_index("gene_id")


def mwu_fisher(frame, resp_col, label):
    """LTR_lnc vs lnc_other on response; Fisher on up fraction."""
    rows = []
    ltr = frame[frame.group3 == "LTR_lnc"][resp_col].dropna()
    oth = frame[frame.group3 == "lnc_other"][resp_col].dropna()
    pcg = frame[frame.group3 == "PCG"][resp_col].dropna()
    if len(ltr) >= 20 and len(oth) >= 20:
        u, p = mannwhitneyu(ltr.values, oth.values, alternative="two-sided")
        rows.append(dict(dataset=label, contrast="LTR_lnc_vs_other",
                         n=len(ltr) + len(oth),
                         median_case=float(np.median(ltr)),
                         median_ctrl=float(np.median(oth)),
                         cliffs=2 * u / (len(ltr) * len(oth)) - 1, pval=p))
        up = resp_col + "_up"
        sub_l = frame[frame.group3 == "LTR_lnc"]
        sub_o = frame[frame.group3 == "lnc_other"]
        a, b = int(sub_l[up].sum()), int(len(sub_l) - sub_l[up].sum())
        c, d = int(sub_o[up].sum()), int(len(sub_o) - sub_o[up].sum())
        orr, pf = fisher_exact([[a, b], [c, d]])
        rows.append(dict(dataset=label, contrast="LTR_lnc_vs_other_UP",
                         n=a + b + c + d,
                         median_case=a / max(a + b, 1),
                         median_ctrl=c / max(c + d, 1),
                         cliffs=float(np.log2(orr)) if orr > 0 else np.nan,
                         pval=pf))
    if len(pcg) >= 20:
        rows.append(dict(dataset=label, contrast="PCG_median_reference",
                         n=len(pcg), median_case=float(np.median(pcg)),
                         median_ctrl=np.nan, cliffs=np.nan, pval=np.nan))
    return rows


def group3_col(frame):
    return "group3"


def cpm_log(mat):
    lib = mat.sum(axis=0)
    return np.log2(mat.div(lib, axis=1) * 1e6 + 1)


def prep_response(counts, cols_treat, cols_ctrl, label, gene_col="gene_id"):
    sub_t = counts[cols_treat]
    sub_c = counts[cols_ctrl]
    keep = (sub_t.mean(axis=1) + sub_c.mean(axis=1)) > 0
    counts = counts[keep].copy()
    lg = cpm_log(counts[cols_treat + cols_ctrl])
    resp = lg[cols_treat].mean(axis=1) - lg[cols_ctrl].mean(axis=1)
    detected = lg[cols_treat + cols_ctrl].mean(axis=1) >= 1.0
    out = pd.DataFrame({"gene_id": counts[gene_col], "response": resp,
                        "detected": detected})
    out["response_up"] = (out.response > 1) & out.detected
    return out


def base_id(s):
    return s.astype(str).str.split(".").str[0]


def gse85461_fastq(rows, med):
    """GSE85461 from our own salmon quants (MB135 iDUX4 + lenti-DUX4 arms)."""
    qd = os.path.join(PD, "dux4", "quants")
    t2g = pd.read_csv(os.path.join(WORKDIR, "data", "gbm_own",
                                   "tx2gene.tsv"), sep="\t")
    t2g = t2g.rename(columns={t2g.columns[0]: "Name",
                              t2g.columns[1]: "gene_id"})
    def gene_tpm(run):
        q = pd.read_csv(os.path.join(qd, run, "quant.sf"), sep="\t")
        q = q.merge(t2g, on="Name", how="inner")
        g = q.groupby("gene_id")["TPM"].sum()
        g.name = run
        return g
    runs_nodox = ["SRR4019004", "SRR4019006", "SRR4019008"]
    runs_dox = ["SRR4019005", "SRR4019007", "SRR4019009"]
    mat = pd.concat([gene_tpm(r) for r in runs_nodox + runs_dox]
                    + [gene_tpm("lenti_Ctrl"), gene_tpm("lenti_DUX4")],
                    axis=1).fillna(0)
    lg = np.log2(mat + 1)
    lab = load_classes().reset_index()
    for label, t, k in [
            ("GSE85461_iDUX4_fastq", runs_dox, runs_nodox),
            ("GSE85461_lentiDUX4_14h", ["lenti_DUX4"], ["lenti_Ctrl"])]:
        resp = lg[t].mean(axis=1) - lg[k].mean(axis=1)
        d = pd.DataFrame({"gene_id": mat.index.values,
                          "response": resp.values,
                          "detected": (lg[t + k].mean(axis=1) >= 1.0).values})
        d["response_up"] = (d.response > 1) & d.detected
        f = d.merge(lab, on="gene_id")
        rows += mwu_fisher(f, "response", label)
        for g in ["LTR_lnc", "lnc_other", "PCG"]:
            v = f[f.group3 == g]["response"].dropna()
            med.append(dict(dataset=label, group=g, n=len(v),
                            median_response=float(np.median(v))))
        tgx = pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t",
                          index_col=0)
        symmap = pd.Series(tgx["symbol"].values,
                           index=tgx["gene_id"].values).reindex(mat.index) \
            .fillna("")
        for g in ["USP29", "TRIM43", "ZSCAN4", "MBD3L2", "RFPL1"]:
            hit = [i for i, s in zip(mat.index, symmap) if s == g][:1]
            if hit:
                print("  sanity", label, g, round(float(resp.loc[hit[0]]), 2))
    return rows


def main():
    with step("30_perturb") as st:
        lab = load_classes()
        rows = []
        med = []

        # ---- 1a. GSE85461 iDUX4 (provided logFC, base ENSG) ---------------
        fc = pd.read_csv(os.path.join(PD, "dux4",
                                      "GSE85461_foldchange_iDUX4.csv.gz"))
        fc["gene_id"] = fc["id"].map(
            pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t",
                        index_col=0)["gene_id"])
        f1 = fc.dropna(subset=["gene_id"]).merge(lab, on="gene_id")
        f1 = f1.assign(response=f1.iDUX4_logFC,
                       response_up=(f1.iDUX4_logFC > 1)
                       & (f1.iDUX4_logCPM > 1))
        rows += mwu_fisher(f1, "response", "GSE85461_iDUX4")

        # ---- 1b. GSE56787 FSHD vs control (counts, base ENSG) -------------
        with gzip.open(os.path.join(PD, "dux4",
                                    "GSE56787_cultured.gene.counts.csv.gz"),
                       "rt") as fh:
            c = pd.read_csv(fh)
        c["gene_id"] = c["gene_id"].map(
            pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t",
                        index_col=0)["gene_id"])
        c = c.dropna(subset=["gene_id"]).drop_duplicates("gene_id") \
            .set_index("gene_id")
        # myotubes only ("tubes"): FSHD vs control, same differentiation state
        t = [x for x in c.columns
             if x.startswith("F") and x.endswith("_tubes")]
        k = [x for x in c.columns
             if x.startswith("C") and x.endswith("_tubes")]
        r = prep_response(c.reset_index(), t, k, "x")
        f2 = r.merge(lab.reset_index(), on="gene_id")
        rows += mwu_fisher(f2, "response", "GSE56787_FSHD_vs_ctrl")

        # ---- 2a. GSE337536 decitabine melanoma (base ENSG) -----------------
        with gzip.open(os.path.join(
                PD, "aza", "GSE337536_Raw_Counts_All_Samples.csv.gz"),
                "rt") as fh:
            a = pd.read_csv(fh)
        a = a.rename(columns={"ID": "ensg"})
        a["gene_id"] = a["ensg"].map(
            pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t",
                        index_col=0)["gene_id"])
        a = a.dropna(subset=["gene_id"]).set_index("gene_id")
        t = [x for x in a.columns if "DEC" in x and "D21" not in x]
        k = [x for x in a.columns if "UNT" in x]
        r = prep_response(a.reset_index(), t, k, "x")
        f3 = r.merge(lab.reset_index(), on="gene_id")
        rows += mwu_fisher(f3, "response", "GSE337536_decitabine")

        # ---- 2b. GSE305426 MPNST DAC/AZA vs DMSO --------------------------
        with gzip.open(os.path.join(
                PD, "aza", "GSE305426_subread_counts.txt.gz"), "rt") as fh:
            m = pd.read_csv(fh, sep="\t", comment="#")
        m = m.rename(columns={"Geneid": "ensg"})
        m["gene_id"] = base_id(m["ensg"]).map(
            pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t",
                        index_col=0)["gene_id"])
        m = m.dropna(subset=["gene_id"]).set_index("gene_id")
        m = m[[x for x in m.columns if not x.startswith(("Chr", "Start",
                                                          "End", "Strand",
                                                          "Length"))]]
        t = [x for x in m.columns if "-DAC-" in x or "-AZA-" in x]
        k = [x for x in m.columns if "-DMSO-" in x]
        r = prep_response(m.reset_index(), t, k, "x")
        f4 = r.merge(lab.reset_index(), on="gene_id")
        rows += mwu_fisher(f4, "response", "GSE305426_DNMTi_vs_DMSO")

        # ---- 1c. GSE85461 from own salmon quants (if chain completed) -----
        if os.path.isdir(os.path.join(PD, "dux4", "quants")) and \
                os.path.exists(os.path.join(PD, "dux4", "quants",
                                            "lenti_DUX4", "quant.sf")):
            rows = gse85461_fastq(rows, med)

        T = pd.DataFrame(rows)
        T.to_csv(safe_write_path(os.path.join(OUT, "perturb_tests.tsv")),
                 sep="\t", index=False)
        for label, f in [("GSE85461_iDUX4", f1), ("GSE56787_FSHD_vs_ctrl",
                                                  f2),
                         ("GSE337536_decitabine", f3),
                         ("GSE305426_DNMTi_vs_DMSO", f4)]:
            for g in ["LTR_lnc", "lnc_other", "PCG"]:
                v = f[f.group3 == g]["response"].dropna()
                med.append(dict(dataset=label, group=g, n=len(v),
                                median_response=float(np.median(v))))
        M = pd.DataFrame(med)
        M.to_csv(safe_write_path(os.path.join(OUT, "response_medians.tsv")),
                 sep="\t", index=False)
        key = T[T.contrast == "LTR_lnc_vs_other"]
        st.note("perturb LTR-vs-other: " + "; ".join(
            "{}: med {} v {}, p={:.2e}".format(r.dataset, r.median_case,
                                               r.median_ctrl, r.pval)
            for r in key.itertuples()))
        st.summary_md("results/perturb/step_summaries/30_perturb.md",
                      "30 perturbation validation (DUX4 + DNMTi)",
                      "## LTR-derived lnc response vs other lnc\n\n"
                      + df_to_md(T, ".4g")
                      + "\n\n## Group median response per dataset\n\n"
                      + df_to_md(M, ".3f"))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    main()
