#!/usr/bin/env python
"""13b_parse_toil.py — Toil matrix × promoter classes + METADATA CROSSWALK
+ SAMPLE ACCOUNTING (per user requirement: track raw-volume2 metadata mapping
and exact per-dataset sample counts used).

Outputs (data/tissue_disease/):
  toil_class_expr.npz   log2(TPM+1) matrix: class genes x samples
  toil_genes.tsv        gene annotation (id/symbol/class/tss_in_ltr/cage_active)
  toil_samples.tsv      sample annotation (cohort/tissue/cancer/condition)
  sample_counts.tsv     ACCOUNTING: per cohort/tissue/cancer n_used (+filters)
  crosswalk_tcga.tsv    volume2 raw-file UUID -> barcode -> cancer -> stage
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

TD = os.path.join(WORKDIR, "data", "tissue_disease")
RES = os.path.join(WORKDIR, "results")


def barcode_type(s):
    if not s.startswith("TCGA"):
        return "other"
    parts = s.split("-")
    code = parts[3][:2] if len(parts) > 3 else ""
    if code.startswith("01") or code.startswith("06"):
        return "tumor"
    if code.startswith("11"):
        return "normal"
    return "other"


def main():
    with step("13b_parse_toil") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                           index_col=0, low_memory=False)
        meta["base_id"] = [g.split(".")[0] for g in meta.index]
        keep = set(meta["base_id"])

        # ---- read matrix
        row_names, chunks, buf = [], [], []
        with gzip.open(os.path.join(TD, "toil_tpm.gz"), "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                gid = f[0].split(".")[0]
                if gid in keep:
                    row_names.append(gid)
                    buf.append(f[1:])
                if len(buf) >= 5000:
                    chunks.append(np.array(buf, dtype=np.float32))
                    buf = []
        if buf:
            chunks.append(np.array(buf, dtype=np.float32))
        M = np.log2(np.vstack(chunks) + 1)
        samples = header[1:]
        st.note(f"matrix: {M.shape[0]} class-genes x {len(samples)} samples")

        # ---- sample annotation
        samp = pd.DataFrame({"sample": samples})
        samp["cohort"] = np.where(samp["sample"].str.upper().str.startswith("GTEX"),
                                  "GTEx", "TCGA")
        satt = pd.read_csv(os.path.join(TD, "gtex_sample_attributes.txt"),
                           sep="\t", skiprows=10, low_memory=False)
        smtsd = [c for c in satt.columns if c.upper().endswith("SMTSD")][0]
        sm = satt[["SAMPID", smtsd]].dropna()
        sm.columns = ["sample", "tissue"]
        samp = samp.merge(sm, on="sample", how="left")
        samp["condition"] = samp["sample"].map(barcode_type)
        samp.loc[samp.cohort == "GTEx", "condition"] = "normal_tissue"

        ss = pd.read_csv(os.path.join(TD, "gdc_sample_sheet.2024-07-09.tsv"),
                         sep="\t")
        ss["sample15"] = ss["Sample ID"].str[:15]
        proj = ss[["sample15", "Project ID"]].drop_duplicates("sample15")
        samp["cancer"] = samp["sample"].str[:15].map(
            proj.set_index("sample15")["Project ID"]).values

        # ---- CROSSWALK: volume2 raw BAM uuid -> barcode -> cancer -> stage
        stage = pd.read_csv(os.path.join(TD, "stage_anno.txt"), sep="\t")
        stage.columns = ["file_name", "sample_type", "project", "stage"]
        stage["file_id"] = stage["file_name"].str.replace(
            r"\.bam.*$", "", regex=True)
        ss2 = ss.rename(columns={"File ID": "file_id", "Sample ID": "barcode",
                                 "Project ID": "project"})
        cw = ss2[["file_id", "barcode", "project"]].merge(
            stage[["file_id", "sample_type", "project", "stage"]],
            on="file_id", how="outer", suffixes=("", "_stage"))
        cw["patient12"] = cw["barcode"].str[:12]
        cw.to_csv(safe_write_path(os.path.join(TD, "crosswalk_tcga.tsv")),
                  sep="\t", index=False)
        st.note(f"crosswalk rows={len(cw)}; with stage="
                f"{int(cw['stage'].notna().sum())}")

        # ---- genes table
        mbi = meta.set_index("base_id")
        base2full = mbi.index.map(lambda b: meta.loc[meta.base_id == b].index[0])
        genes = pd.DataFrame(index=row_names)
        genes["gene_id"] = [meta.index[meta.base_id == b][0] for b in row_names]
        genes["symbol"] = mbi["gene_name"].reindex(row_names).values
        genes["cls"] = mbi["cls"].reindex(row_names).values
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                         index_col=0, low_memory=False)[["tss_in_ltr"]]
        act = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                          index_col=0, low_memory=False)[["active"]]
        te["base_id"] = [g.split(".")[0] for g in te.index]
        act["base_id"] = [g.split(".")[0] for g in act.index]
        genes["tss_in_ltr"] = te.set_index("base_id")["tss_in_ltr"].reindex(
            row_names).values
        genes["cage_active"] = act.set_index("base_id")["active"].reindex(
            row_names).values

        np.savez_compressed(safe_write_path(
            os.path.join(TD, "toil_class_expr.npz")), M=M)
        genes.to_csv(safe_write_path(os.path.join(TD, "toil_genes.tsv")), sep="\t")
        samp.to_csv(safe_write_path(os.path.join(TD, "toil_samples.tsv")),
                    sep="\t", index=False)

        # ---- SAMPLE ACCOUNTING
        acct = []
        acct.append(dict(dataset="Toil matrix", field="all_samples",
                         n=len(samp), note="columns in matrix"))
        acct.append(dict(dataset="Toil matrix", field="class-genes matched",
                         n=len(genes), note="of 63,086 GENCODE genes"))
        g = samp[samp.cohort == "GTEx"]
        acct.append(dict(dataset="GTEx", field="samples", n=len(g),
                         note="with tissue annotation: "
                              f"{int(g.tissue.notna().sum())}"))
        for t, sub in g.groupby("tissue"):
            if len(sub) >= 20:
                acct.append(dict(dataset="GTEx", field=t, n=len(sub), note="tissue"))
        t_ = samp[samp.cohort == "TCGA"]
        acct.append(dict(dataset="TCGA", field="tumor samples",
                         n=int((t_.condition == "tumor").sum()), note="barcode 01/06"))
        acct.append(dict(dataset="TCGA", field="normal samples",
                         n=int((t_.condition == "normal").sum()), note="barcode 11"))
        for c, sub in t_[t_.condition == "tumor"].groupby("cancer"):
            nn = int(((samp.cancer == c) & (samp.condition == "normal")).sum())
            if len(sub) >= 20:
                acct.append(dict(dataset="TCGA", field=f"{c} tumor",
                                 n=len(sub),
                                 note=f"paired normals={nn}"))
        ac = pd.DataFrame(acct)
        ac.to_csv(safe_write_path(os.path.join(TD, "sample_counts.tsv")),
                  sep="\t", index=False)
        st.summary_md("results/tissue_disease/step_summaries/13b_sample_counts.md",
                      "13b sample accounting",
                      "## Samples used per dataset/tissue/cancer\n\n"
                      + df_to_md(ac, ".0f", max_rows=80))
        st.note(f"accounting rows={len(ac)}; GTEx tissues(n>=20)="
                f"{sum(1 for r in acct if r['dataset']=='GTEx' and r['field']!='samples')};")


if __name__ == "__main__":
    main()
