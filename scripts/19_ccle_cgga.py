#!/usr/bin/env python
"""19_ccle_cgga.py — CCLE cell-line & CGGA glioma extension (actual formats).

CCLE (DepMap new portal format):
  rows = cell lines (SequencingID, ModelConditionID, ModelID, ...),
  expr columns = "SYMBOL (EntrezID)"; lineage from OmicsProfiles.csv.
CGGA (mRNAseq_325 + _693):
  rows = Gene_Name symbols, cols = CGGA_ID, values = RSEM TPM;
  clinical: Grade (WHO II/III/IV), IDH_mutation_status, Histology.

Analyses:
  CCLE: class activity by lineage; DUX4 expression vs LTR-derived lncRNA
        activation (route-model test: Spearman + decile contrast).
  CGGA: class activity & LTR derepression by grade and IDH status
        (validates the TCGA-GBM LTR signal at n≈1,000).
Outputs -> results/ccle_cgga/tables/*.tsv + step summary
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

TD = os.path.join(WORKDIR, "data", "ccle_cgga")
OUT = os.path.join(WORKDIR, "results", "ccle_cgga", "tables")
RES = os.path.join(WORKDIR, "results")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]


def gene_maps():
    meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                       low_memory=False)
    sym2cls = dict(zip(meta.gene_name.astype(str), meta.cls))
    te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                     index_col=0, low_memory=False)[["tss_in_ltr"]]
    te = te.join(meta.set_index("gene_id")["gene_name"])
    sym2ltr = {str(s): int(v) for s, v in zip(te.gene_name, te.tss_in_ltr)
               if pd.notna(s)}
    return sym2cls, sym2ltr


def log2tpm(v):
    """DataFrame/Series of TPM (or already log2) -> log2(x+1) if needed."""
    if isinstance(v, pd.DataFrame):
        v = v.apply(pd.to_numeric, errors="coerce")
        return v if float(v.max().max()) < 25 else np.log2(v + 1)
    v = pd.to_numeric(pd.Series(v), errors="coerce")
    return v if float(v.max()) < 25 else np.log2(v + 1)


# --------------------------------------------------------------- CCLE -------
def ccle(sym2cls, sym2ltr):
    p_expr = os.path.join(TD, "ccle", "OmicsExpressionTPMLogp1HumanAllGenes.csv")
    p_prof = os.path.join(TD, "ccle", "OmicsProfiles.csv")
    if not (os.path.exists(p_expr) and os.path.exists(p_prof)):
        return "skip: CCLE matrices absent"
    d = pd.read_csv(p_expr, index_col=0, low_memory=False)
    meta_cols = [c for c in ["SequencingID", "ModelConditionID", "ModelID",
                             "IsDefaultEntryForMC", "IsDefaultEntryForModel"]
                 if c in d.columns]
    d = d.drop_duplicates(subset=[c for c in meta_cols if c == "ModelID"]
                          or [d.columns[0]])
    prof = pd.read_csv(p_prof)
    prof = prof[(prof.DataType == "rna") & (prof.IsDefaultEntryForMC.astype(str).isin(["Yes","True"]))]
    m2lineage = prof.drop_duplicates("ModelID").set_index("ModelID")["Lineage"]
    d["Lineage"] = d["ModelID"].map(m2lineage) if "ModelID" in d else "NA"

    expr_cols = [c for c in d.columns if "(" in c]
    syms = [c.split(" (")[0] for c in expr_cols]
    cls_of = {c: sym2cls.get(s) for c, s in zip(expr_cols, syms)}
    ltr_of = {c: sym2ltr.get(s, 0) for c, s in zip(expr_cols, syms)}
    dux4_col = next((c for c in expr_cols if c.startswith("DUX4 (")), None)

    rows = []
    for lineage, sub in d.groupby("Lineage"):
        if len(sub) < 15:
            continue
        for c in CLASSES + ["PCG"]:
            cols = [k for k, v in cls_of.items() if v == c]
            vals = log2tpm(sub[cols])
            rows.append(dict(lineage=lineage, n_lines=len(sub), cls=c,
                             frac_expr=float((vals >= 1).any(axis=1).mean()),
                             med_active=float(np.nanmedian(
                                 vals[(vals >= 1).any(axis=1)].mean(axis=1)))
                             if vals.size else np.nan))
        lcols = [k for k in expr_cols if cls_of.get(k) in CLASSES and ltr_of[k] == 1]
        if lcols:
            vals = log2tpm(sub[lcols])
            rows.append(dict(lineage=lineage, n_lines=len(sub),
                             cls="LTR_derived_lnc",
                             frac_expr=float((vals >= 1).any(axis=1).mean()),
                             med_active=float(np.nanmedian(
                                 vals[(vals >= 1).any(axis=1)].mean(axis=1)))))
    t1 = pd.DataFrame(rows)
    t1.to_csv(safe_write_path(os.path.join(OUT, "ccle_class_by_lineage.tsv")),
              sep="\t", index=False)

    # DUX4 vs LTR-lncRNA activation
    res = None
    if dux4_col is not None:
        lcols = [k for k in expr_cols if cls_of.get(k) in CLASSES and ltr_of[k] == 1]
        ocols = [k for k in expr_cols if cls_of.get(k) in CLASSES and ltr_of[k] == 0]
        ltr_act = log2tpm(d[lcols]).mean(axis=1)
        oth_act = log2tpm(d[ocols]).mean(axis=1)
        rel = ltr_act - oth_act
        dux = log2tpm(d[dux4_col])
        print(f"[diag] DUX4 col={dux4_col}; lncLTR n={len(lcols)}; "
              f"DUX4>1 in {int((dux > 1).sum())}/{len(dux)} lines; "
              f"DUX4 max={float(dux.max()):.1f}")
        ok = dux.notna() & rel.notna()
        rho, p = spearmanr(dux[ok], rel[ok])
        hi = ok & (dux >= dux[ok].quantile(0.9))
        lo = ok & (dux <= dux[ok].quantile(0.1))
        u, pu = mannwhitneyu(rel[hi], rel[lo])
        res = dict(rho=rho, p_spearman=p, hi_minus_lo=float(rel[hi].mean() - rel[lo].mean()),
                   p_mwu=pu, n=int(ok.sum()))
        pd.DataFrame([res]).to_csv(safe_write_path(os.path.join(
            OUT, "ccle_dux4_vs_ltr.tsv")), sep="\t", index=False)
    msg = (f"CCLE: {t1.lineage.nunique()} lineages, {len(d)} lines; "
           f"DUX4~LTR-lnc rel-activation rho={res['rho']:.2f} "
           f"(p={res['p_spearman']:.1e}, hi-lo={res['hi_minus_lo']:+.2f}, "
           f"MWU p={res['p_mwu']:.1e}, n={res['n']})" if res else
           f"CCLE: {len(d)} lines (DUX4 column not found)")
    return msg, res


# --------------------------------------------------------------- CGGA -------
def cgga(sym2cls, sym2ltr):
    frames, clins = [], []
    for cohort in ("CGGA_mRNAseq_325", "CGGA_mRNAseq_693"):
        p = os.path.join(TD, "cgga", cohort, f"{cohort.replace('CGGA_', 'CGGA.').replace('CGGA.', 'CGGA.mRNAseq_')}.RSEM-genes.20200506.txt")
        p = os.path.join(TD, "cgga", cohort,
                         cohort.replace("CGGA_mRNAseq", "CGGA.mRNAseq")
                         + ".RSEM-genes.20200506.txt")
        pc = p.replace(".RSEM-genes.20200506.txt", "_clinical.20200506.txt")
        if not (os.path.exists(p) and os.path.exists(pc)):
            continue
        m = pd.read_csv(p, sep="\t", index_col=0)
        frames.append(m)
        clins.append(pd.read_csv(pc, sep="\t"))
    if not frames:
        return "skip: CGGA matrices absent", None
    d = pd.concat(frames, axis=1)
    d = d.loc[~d.index.duplicated()]
    clin = pd.concat(clins, ignore_index=True)
    clin = clin.drop_duplicates(subset=["CGGA_ID"]).set_index("CGGA_ID")
    common = d.columns.intersection(clin.index)
    d = d[common]
    clin = clin.loc[common]

    d.index = d.index.astype(str)
    cls = np.array([sym2cls.get(g) for g in d.index])
    ltr = np.array([sym2ltr.get(g, 0) for g in d.index])
    d = d.apply(pd.to_numeric, errors="coerce")
    dl = log2tpm(d)

    rows = []
    for grp_name, col in [("Grade", "Grade"), ("IDH", "IDH_mutation_status")]:
        if col not in clin.columns:
            continue
        for grp, idx in clin.groupby(col).groups.items():
            idx = [i for i in idx if i in dl.columns]
            if len(idx) < 15:
                continue
            block = dl[idx]
            for c in CLASSES + ["PCG"]:
                sub = block[cls == c]
                if sub.shape[0] < 10:
                    continue
                rows.append(dict(group=grp_name, level=str(grp), n=len(idx),
                                 cls=c, n_genes=sub.shape[0],
                                 frac_expr=float((sub >= 1).any(axis=1).mean()),
                                 med=float(np.nanmedian(sub.mean(axis=1)))))
            lm = (cls.isin(CLASSES) if hasattr(cls, "isin") else np.isin(cls, CLASSES)) & (ltr == 1)
            om = np.isin(cls, CLASSES) & (ltr == 0)
            if lm.sum() >= 10 and om.sum() >= 50:
                lmed = float(np.nanmedian(block[lm].mean(axis=1)))
                omed = float(np.nanmedian(block[om].mean(axis=1)))
                rows.append(dict(group=grp_name, level=str(grp), n=len(idx),
                                 cls="LTR_derived_lnc", n_genes=int(lm.sum()),
                                 frac_expr=float((block[lm] >= 1).any(axis=1).mean()),
                                 med=lmed, ltr_minus_other=lmed - omed))
    t = pd.DataFrame(rows)
    t.to_csv(safe_write_path(os.path.join(OUT, "cgga_class_by_group.tsv")),
             sep="\t", index=False)
    ltr_rows = t[t.cls == "LTR_derived_lnc"]
    msg = (f"CGGA: {len(common)} samples; "
           + "; ".join(f"{r.group} {r.level}: LTR-Δ={r.ltr_minus_other:+.2f}"
                       for r in ltr_rows.itertuples()))
    return msg, t


def main():
    with step("19_ccle_cgga") as st:
        os.makedirs(OUT, exist_ok=True)
        sym2cls, sym2ltr = gene_maps()
        msg1, ccle_res = ccle(sym2cls, sym2ltr)
        st.note(msg1)
        msg2, cgga_t = cgga(sym2cls, sym2ltr)
        st.note(msg2)
        body = f"- {msg1}\n- {msg2}"
        if cgga_t is not None and len(cgga_t):
            body += "\n\n" + df_to_md(cgga_t.head(40), ".3f", max_rows=40)
        ccle_t = os.path.join(OUT, "ccle_class_by_lineage.tsv")
        if os.path.exists(ccle_t):
            body += ("\n\n## CCLE lineage table (head)\n\n" +
                     df_to_md(pd.read_csv(ccle_t, sep="\t").head(30), ".3f"))
        st.summary_md("results/ccle_cgga/step_summaries/19_ccle_cgga.md",
                      "19 CCLE/CGGA extension", body)


if __name__ == "__main__":
    main()
