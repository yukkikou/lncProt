# 25 Sample accounting (Supplementary Table)

_Generated 2026-10-08 16:40:02 by 25_sample_accounting_

All counts derived live from primary project tables (constants flagged in note/source).

| dataset | unit | item | n | note | source |
|---|---|---|---|---|---|
| GENCODE v46 (project master) | genes | lncRNA classified (5 birth classes) | 19258 | divergent 2692; antisense 2778; intronic 3403; sense_overlapping 350; intergenic 10035 | results/gene_master_table.tsv |
| GENCODE v46 (project master) | genes | PCG baseline | 20065 | CAGE-refined TSS 91% | results/gene_master_table.tsv |
| GENCODE v46 (project master) | genes | pseudogenes (incl. transcribed) | 15210 | comparator group | results/gene_master_table.tsv |
| GENCODE v46 (project master) | genes | genes matched in Toil matrix | 48275 | of 63,086 GENCODE genes | data/tissue_disease/toil_genes.tsv |
| FANTOM5 CAGE (hg38 fair+new) | samples | TPM matrix samples | 1829 | matrix columns minus 7 annotation columns | data/fantom/hg38_cage_tpm_ann.osc.txt.gz |
| Roadmap ChromHMM (hg38-lifted) | cell types | segmentations used | 5 | E003/E017/E067/E071/E088 | data/chromhmm_ov_*.tsv |
| HOCOMOCO v11 + JASPAR2024 (FIMO) | TF motifs | TFs with promoter hits (presence matrix) | 1131 | 1,176 motifs scanned | results/motif_presence.tsv |
| ReMap 2022 | TFs | TFs with promoter binding | 1234 | 38.4M binding records | results/tf_binding_enrichment.tsv |
| Candidate TFs (dual evidence) | TFs | motif+binding agreement | 91 | class-specific rows | results/candidate_TFs.tsv |
| Xena Toil (TCGA+GTEx unified) | samples | matrix samples | 19131 | 48,275 class genes matched | data/tissue_disease/toil_samples.tsv |
| GTEx (within Toil) | samples | samples with tissue annotation, n>=20 tissues | 7741 | 49 tissues | data/tissue_disease/toil_samples.tsv |
| TCGA (within Toil) | samples | tumour samples (barcode 01/06) | 9578 | 33 cancer types in matrix | data/tissue_disease/toil_samples.tsv |
| TCGA (within Toil) | samples | normal samples (barcode 11) | 727 | usable contrasts in cancers with >=5 normals | data/tissue_disease/toil_samples.tsv |
| TCGA (within Toil) | samples | tumour+normal in >=5-normal cancers | 7348 | 18 cancers (Fig8a set) | data/tissue_disease/toil_samples.tsv |
| CGGA glioma mRNA-seq | samples | Grade levels | 1013 | WHO II 291; WHO III 334; WHO IV 388 | results/ccle_cgga/tables/cgga_class_by_group.tsv |
| CGGA glioma mRNA-seq | samples | IDH levels | 966 | Mutant 531; Wildtype 435 | results/ccle_cgga/tables/cgga_class_by_group.tsv |
| CGGA glioma mRNA-seq | samples | total profiled | 1018 | packages mRNAseq_325 + mRNAseq_693 (constant: module-19 log) | logs/PROGRESS.md |
| EPIC methylation (IDH-A) | samples | paired with RNA-seq | 28 | 765,084 probes; 29 samples processed (constant) | results/ccle_cgga/tables/methylation_vs_ltr.tsv |
| scRNA-seq (IDH-A) | cells | typed cells tabulated | 727 | Astrocyte 66; Myeloid 101; Neoplastic 417; Oligodendrocyte 143; of 765 total cells (38 LowQC/T/B not tabulated); 18 patients (constant: module-23) | results/scrna/class_by_celltype.tsv |
| CCLE (DepMap) | cell lines | lines with class detection | 1666 | 24 lineages; of 1,719 lines in matrix (constant: module-19) | results/ccle_cgga/tables/ccle_class_by_lineage.tsv |
| Medulloblastoma cohort (public) | libraries | salmon-quantified libraries | 206 | heterogeneous identifiers | data/mb_own/quants/*/quant.sf |
| Medulloblastoma cohort (public) | samples | samples after library aggregation | 175 | stem-core aggregation; no clinical annotation available | results/cns_cohorts/mb_per_sample.tsv |
