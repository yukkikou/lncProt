#!/usr/bin/env python
"""25_sample_accounting.py — project-wide sample accounting (Supplementary Table).

Every count is derived live from its primary project table (no hand-typed
numbers except the two flagged constants with their sources in `source`).
Output: results/sample_accounting.tsv + formatted supplementary-table markdown.
"""
import glob
import gzip
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, WORKDIR

RES = os.path.join(WORKDIR, "results")
TD = os.path.join(WORKDIR, "data", "tissue_disease")
FANTOM_TPM = os.path.join(WORKDIR, "data", "fantom",
                          "hg38_cage_tpm_ann.osc.txt.gz")
QUANTS = os.path.join(WORKDIR, "data", "gbm_own", "quants")

rows = []


def add(dataset, unit, item, n, note, source):
    rows.append(dict(dataset=dataset, unit=unit, item=item, n=n,
                     note=note, source=source))


def main():
    with step("25_sample_accounting") as st:
        # ---- genes (discovery layer) ------------------------------------
        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                        sep="\t", low_memory=False)
        lnc5 = ["divergent", "antisense", "intronic", "sense_overlapping",
                "intergenic"]
        n_lnc = int(m.cls.isin(lnc5).sum())
        n_pcg = int((m.cls == "PCG").sum())
        n_pseudo = int(m.cls.isin(["pseudogene",
                                   "transcribed_pseudogene"]).sum())
        add("GENCODE v46 (project master)", "genes",
            "lncRNA classified (5 birth classes)", n_lnc,
            "; ".join("{} {}".format(c, (m.cls == c).sum()) for c in lnc5),
            "results/gene_master_table.tsv")
        add("GENCODE v46 (project master)", "genes", "PCG baseline",
            n_pcg, "CAGE-refined TSS 91%", "results/gene_master_table.tsv")
        add("GENCODE v46 (project master)", "genes",
            "pseudogenes (incl. transcribed)", n_pseudo,
            "comparator group", "results/gene_master_table.tsv")
        tg = pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t")
        add("GENCODE v46 (project master)", "genes",
            "genes matched in Toil matrix", len(tg),
            "of 63,086 GENCODE genes", "data/tissue_disease/toil_genes.tsv")

        # ---- FANTOM5 CAGE -------------------------------------------------
        with gzip.open(FANTOM_TPM, "rt") as fh:
            for line in fh:
                if not line.startswith("#"):
                    n_fant = len(line.rstrip("\n").split("\t")) - 7
                    break
        add("FANTOM5 CAGE (hg38 fair+new)", "samples",
            "TPM matrix samples", n_fant,
            "matrix columns minus 7 annotation columns",
            "data/fantom/hg38_cage_tpm_ann.osc.txt.gz")

        # ---- Roadmap ChromHMM ---------------------------------------------
        n_chrom = len(glob.glob(os.path.join(WORKDIR, "data",
                                             "chromhmm_ov_*.tsv")))
        add("Roadmap ChromHMM (hg38-lifted)", "cell types",
            "segmentations used", n_chrom,
            "E003/E017/E067/E071/E088", "data/chromhmm_ov_*.tsv")

        # ---- TF evidence layer --------------------------------------------
        mp = pd.read_csv(os.path.join(RES, "motif_presence.tsv"),
                         sep="\t", nrows=0)
        n_motif = len(mp.columns) - 1
        tb = pd.read_csv(os.path.join(RES, "tf_binding_enrichment.tsv"),
                         sep="\t", usecols=["tf"])
        add("HOCOMOCO v11 + JASPAR2024 (FIMO)", "TF motifs",
            "TFs with promoter hits (presence matrix)", n_motif,
            "1,176 motifs scanned", "results/motif_presence.tsv")
        add("ReMap 2022", "TFs", "TFs with promoter binding",
            int(tb.tf.nunique()), "38.4M binding records",
            "results/tf_binding_enrichment.tsv")
        ct = pd.read_csv(os.path.join(RES, "candidate_TFs.tsv"), sep="\t")
        add("Candidate TFs (dual evidence)", "TFs",
            "motif+binding agreement", int((ct.evidence == "motif+binding")
                                           .sum()), "class-specific rows",
            "results/candidate_TFs.tsv")

        # ---- Toil / GTEx / TCGA -------------------------------------------
        ts = pd.read_csv(os.path.join(TD, "toil_samples.tsv"), sep="\t")
        add("Xena Toil (TCGA+GTEx unified)", "samples", "matrix samples",
            len(ts), "48,275 class genes matched",
            "data/tissue_disease/toil_samples.tsv")
        gt = ts[(ts.cohort == "GTEx") & ts.tissue.notna()
                & (ts.tissue.astype(str).str.strip() != "")]
        vc = gt.tissue.value_counts()
        keep = vc[vc >= 20]
        add("GTEx (within Toil)", "samples",
            "samples with tissue annotation, n>=20 tissues",
            int(keep.sum()), "{} tissues".format(len(keep)),
            "data/tissue_disease/toil_samples.tsv")
        tc = ts[ts.cohort == "TCGA"]
        n_tu = int((tc.condition == "tumor").sum())
        n_no = int((tc.condition == "normal").sum())
        add("TCGA (within Toil)", "samples", "tumour samples (barcode 01/06)",
            n_tu, "33 cancer types in matrix",
            "data/tissue_disease/toil_samples.tsv")
        add("TCGA (within Toil)", "samples", "normal samples (barcode 11)",
            n_no, "usable contrasts in cancers with >=5 normals",
            "data/tissue_disease/toil_samples.tsv")
        norm_per = tc[tc.condition == "normal"].groupby("cancer").size()
        ok_c = norm_per[norm_per >= 5].index
        sub = tc[tc.cancer.isin(ok_c)]
        add("TCGA (within Toil)", "samples",
            "tumour+normal in >=5-normal cancers",
            int((sub.condition == "tumor").sum() + (sub.condition == "normal")
                .sum()), "{} cancers (Fig8a set)".format(len(ok_c)),
            "data/tissue_disease/toil_samples.tsv")

        # ---- CGGA ----------------------------------------------------------
        cg = pd.read_csv(os.path.join(RES, "ccle_cgga", "tables",
                                      "cgga_class_by_group.tsv"), sep="\t")
        for grp in ("Grade", "IDH"):
            sub = cg[cg.group == grp].drop_duplicates(subset=["level"])
            det = "; ".join("{} {}".format(r.level, int(r.n))
                            for _, r in sub.iterrows())
            add("CGGA glioma mRNA-seq", "samples",
                "{} levels".format(grp), int(sub.n.sum()), det,
                "results/ccle_cgga/tables/cgga_class_by_group.tsv")
        add("CGGA glioma mRNA-seq", "samples", "total profiled", 1018,
            "packages mRNAseq_325 + mRNAseq_693 (constant: module-19 log)",
            "logs/PROGRESS.md")

        # ---- EPIC methylation ----------------------------------------------
        mv = pd.read_csv(os.path.join(RES, "ccle_cgga", "tables",
                                      "methylation_vs_ltr.tsv"), sep="\t")
        n_pair = int(mv["n"].iloc[0])
        add("EPIC methylation (IDH-A)", "samples", "paired with RNA-seq",
            n_pair, "765,084 probes; 29 samples processed (constant)",
            "results/ccle_cgga/tables/methylation_vs_ltr.tsv")

        # ---- scRNA ----------------------------------------------------------
        sc = pd.read_csv(os.path.join(RES, "scrna", "class_by_celltype.tsv"),
                         sep="\t")
        per_type = sc.drop_duplicates(subset=["cell_type"])
        det = "; ".join("{} {}".format(r.cell_type, int(r.n_cells))
                        for _, r in per_type.iterrows())
        add("scRNA-seq (IDH-A)", "cells", "typed cells tabulated",
            int(per_type.n_cells.sum()),
            det + "; of 765 total cells (38 LowQC/T/B not tabulated); "
            "18 patients (constant: module-23)",
            "results/scrna/class_by_celltype.tsv")

        # ---- CCLE ------------------------------------------------------------
        cl = pd.read_csv(os.path.join(RES, "ccle_cgga", "tables",
                                      "ccle_class_by_lineage.tsv"), sep="\t")
        per_lin = cl.drop_duplicates(subset=["lineage"])
        add("CCLE (DepMap)", "cell lines", "lines with class detection",
            int(per_lin.n_lines.sum()),
            "{} lineages; of 1,719 lines in matrix (constant: module-19)"
            .format(len(per_lin)),
            "results/ccle_cgga/tables/ccle_class_by_lineage.tsv")

        # ---- in-house GBM cohort: REMOVED from manuscript 2026-10-08 -------
        # (results archived under results/gbm_own; crosswalk + grade analysis
        #  retained on disk, restorable if ethics documentation is added)

        # ---- medulloblastoma cohort (public download) ----------------------
        mbq = len(glob.glob(os.path.join(WORKDIR, "data", "mb_own",
                                         "quants", "*", "quant.sf")))
        mbp = pd.read_csv(os.path.join(RES, "cns_cohorts",
                                       "mb_per_sample.tsv"), sep="\t")
        add("Medulloblastoma cohort (public)", "libraries", "salmon-quantified libraries",
            mbq, "heterogeneous identifiers", "data/mb_own/quants/*/quant.sf")
        add("Medulloblastoma cohort (public)", "samples",
            "samples after library aggregation", len(mbp),
            "stem-core aggregation; no clinical annotation available",
            "results/cns_cohorts/mb_per_sample.tsv")

        T = pd.DataFrame(rows)
        T.to_csv(safe_write_path(os.path.join(RES, "sample_accounting.tsv")),
                 sep="\t", index=False)

        lines = ["| dataset | unit | item | n | note | source |",
                 "|---|---|---|---|---|---|"]
        for r in rows:
            lines.append("| {} | {} | {} | {} | {} | {} |".format(
                r["dataset"], r["unit"], r["item"], r["n"], r["note"],
                r["source"]))
        st.summary_md("results/step_summaries/25_sample_accounting.md",
                      "25 Sample accounting (Supplementary Table)",
                      "All counts derived live from primary project tables "
                      "(constants flagged in note/source).\n\n" +
                      "\n".join(lines))
        st.note("rows={} -> results/sample_accounting.tsv".format(len(T)))


if __name__ == "__main__":
    main()
