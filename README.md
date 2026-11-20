# lncProt pipeline

Code for "The genomic element supplying the promoter sets promoter architecture but not transcriptional regulation in human lncRNAs".

The pipeline classifies GENCODE v46 lncRNAs by the genomic element that supplies the promoter, scores promoters for sequence, chromatin, conservation, cap activity, factor binding and enhancer RNA, and produces the figures and tables in the manuscript.

## Layout

```
clean_scripts/
├── config.example.sh     copy to config.sh and edit the paths
├── scripts/              analysis modules 01-90, run in numeric order
├── R/                    figure layer (00_common.R, 01-07)
└── figures/              manuscript figure builders (chapter plates, panel redraws, legends)
```

## Requirements

- Python 3.10 with numpy, pandas, scipy, matplotlib, scikit-learn, statsmodels, pyarrow
- R 4.4 with ggplot2, dplyr, patchwork, cowplot, ComplexHeatmap, data.table
- Command-line tools: bedtools 2.30, MEME/FIMO 5.4.1, UCSC liftOver, ghostscript, poppler-utils
- Cluster execution uses SLURM; `scripts/pylib/progress.py` provides the stage wrapper used by every module

Reference genomes, annotations and database mirrors are described in the manuscript's Data Availability section. Public accessions are listed there; no data are redistributed in this repository.

## Setup

```bash
cp config.example.sh config.sh
# edit WORKDIR, REF, DB and the tool paths in config.sh
source config.sh
```

Every module reads paths from `config.sh` through `scripts/pylib/progress.py`. The figure builders in `figures/` carry a `ROOT` constant at the top of the file; set it to the repository path before running them.

## Run order

1. `scripts/01_download.sh` to `scripts/02_parse_gencode.py` — annotation and gene sets
2. `scripts/03a_chain_liftover.py`, `03b_chromhmm.sh`, `04_classify_loci.py` — coordinate classes
3. `scripts/05_seq_features.py`, `06_te_promoters.py`, `07a_fimo_prepare.sh`, `07b_remap_binding.sh`, `08_enrichment.py` — promoter features
4. `scripts/09_cage_activity.py` to `12e_partner_pcgs.py` — activity, conservation, statistics
5. `scripts/13a_download_toil.sh` to `14_td_figures.py` — GTEx and TCGA
6. `scripts/15_final_gaps.py` to `30_perturb.py` — cross-species, methylation, perturbation
7. `scripts/31_tf_contexts.py` to `90_p0_reanalysis.py` — reviewer-driven analyses
8. `scripts/16_submission_figures.py`, `scripts/24_ml_extensions.py` — main figures
9. `R/run_all_plots.R` — R figure layer
10. `figures/` — chapter composites, panel redraws and legends

## Figures

- `figures/extract_panels.py` splits Fig. 5 and Fig. 7 into single panels.
- `figures/redraw_panels.py` redraws the perturbation panel and the three-panel integration figure.
- `figures/build_composites.py` and `figures/build_chapters6.py` assemble the six chapter figures with one uppercase letter per panel.
- `figures/make_legends6.py` builds the chapter legends from the source legends.
- `figures/fig1_v8.py` and `figures/make_fig1B_species.py` build the Figure 1 panels.

## Notes

- The pipeline is restartable. Each module writes a step summary under `results/step_summaries/` and skips completed stages.
- Statistical tests are registered in `results/review_p1/P1_6_test_registry.tsv`.
- The merged annotation is `gene_master_table_merged.tsv`, 39 columns, 54,533 promoters.

## Contact

Ence Yang (yangence@pku.edu.cn), Peking University Health Science Center.
