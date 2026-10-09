# lncProt — Birth-class analysis pipeline

Code for: *Divergent borrowing and endogenous retroviral exaptation define two main routes to transcriptional competence in human lncRNA genes* (Hu & Yang).

## What this does
Coordinate-defined five-way classification of all GENCODE v46 lncRNA loci (divergent / antisense / intronic / sense-overlapping / intergenic, plus pseudogene comparator), uniform promoter scoring (sequence, core elements, TE content, phastCons + 241-mammal phyloP, ChromHMM, CAGE activity, TF motifs/binding), tissue/disease extension (GTEx 49 tissues, TCGA 33 types, CGGA, CCLE, scRNA, public medulloblastoma), perturbation validation (DUX4 overexpression, DNMT inhibitors) and machine-learning analyses.

## Requirements
- Linux + SLURM; ~24 GB RAM for the Toil-matrix modules (28C, 31T)
- Python 3.10 (scipy 1.15, pandas 2.3, scikit-learn 1.7, numpy; no installs needed beyond a conda env with these)
- Singularity/Apptainer containers (read-only binds): salmon 1.10.2, MEME/FIMO 5.4.1, bedtools 2.30.0, pyBigWig, gffread 0.12.7, R/Seurat

## Quick start
```bash
bash slurm/run_all.sh          # marker-file restartable: re-run skips completed stages
python scripts/28_evo_exts.py A   # any single module can be invoked directly
```

## Module map (01–31 → paper sections)
| Modules | Content | Key outputs |
|---|---|---|
| 01, 13a | reference/data acquisition | data/fantom, Toil matrix |
| 02–04b | parse GENCODE, classify loci | results/gene_master_table.tsv |
| 05–06 | promoter sequence + TE features | seq_features, te_promoter_stats |
| 07a–08 | FIMO motifs, ReMap binding, enrichment | candidate_TFs.tsv |
| 09–10, 27 | CAGE activity, phastCons/ChromHMM, phyloP+acceleration | results/phylop/ |
| 12a–12e | publication statistics + robustness | pub_stats, robustness_* |
| 13b–13e, 14, 15 | GTEx/TCGA tissue & disease | results/tissue_disease/ |
| 19–23, s7–s9 | CCLE/CGGA/EPIC/scRNA validation | results/ccle_cgga/, scrna/ |
| 22, 26, s10–s11 | own GBM (archived) + public MB cohorts | results/gbm_own/, cns_cohorts/ |
| 24, 28 | ML extensions, phyloP ablation/tiers | results/ml/, results/evo_ext/ |
| 25 | sample accounting | results/sample_accounting.tsv |
| 29, 30, 31 | promoter shape/age/η², perturbation, TF contexts | promoter_ext/, perturb/, tf_contexts/ |
| 11, 16 | report + submission figures | REPORT.md, figures/sub/ |

## Reproducibility guarantees
- Fixed seeds (42) everywhere; ML folds deterministic; **hard gates**: module 28A aborts unless baseline macro-F1 reproduces 0.336; module 28C aborts unless per-gene TCGA DE matches the 13d shards to <1e-9; fastq MD5s checked before quantification.
- Full run ledger auto-written to `logs/PROGRESS.md`; per-module summaries under `results/**/step_summaries/`.

## Data sources (all public)
GENCODE v46; FANTOM5 hg38 (fair+new); GTEx dbGaP phs000424.v8.p2; TCGA via UCSC Xena Toil; ReMap 2022; Roadmap ChromHMM (hg19, lifted in-pipeline); HOCOMOCO v11; JASPAR 2024; phastCons 470-way + Cactus 241-mammal phyloP (UCSC); CGGA; CCLE/DepMap; GEO GSE85461/GSE56787/GSE337536/GSE305426; GSA-Human HRA011170.

## Package contents
`scripts/` (31 modules + pylib), `slurm/` (s15 sbatch drivers), `config.sh` (paths), `results/step_summaries/`, `logs/PROGRESS.md`, `REFERENCES.md`, `SUPPLEMENTARY_INDEX.md`.
