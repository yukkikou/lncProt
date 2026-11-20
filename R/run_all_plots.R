# =============================================================================
# run_all_plots.R  --  run every figure script in order
# -----------------------------------------------------------------------------
# Usage (from the lncProt project root):
#   Rscript R/run_all_plots.R
# or interactively:
#   source("R/run_all_plots.R")
#
# Figures are written to figures/R/ as vector PDF and 600-dpi PNG.
# Required packages: tidyverse, patchwork, scales.
# =============================================================================

PROJ <- Sys.getenv("LNCPROT_DIR", unset = normalizePath("."))
Sys.setenv(LNCPROT_DIR = PROJ)
if (!dir.exists(file.path(PROJ, "results"))) {
  stop("results/ not found under ", PROJ,
       ". Run from the lncProt project root or set LNCPROT_DIR.")
}

scripts <- c(
  "00_common.R",                 # shared palettes, theme, helpers
  "01_classes_promoters.R",      # figures 01-05
  "02_te_conservation_tf.R",     # figures 06-13
  "03_activity_constraint.R",    # figures 14-19
  "04_tissue_cancer.R",          # figures 20-26
  "05_evolution_parity.R",       # figures 27-33
  "06_ml_perturbation.R",        # figures 34-37
  "07_composite_figures.R"       # composites C1-C6 (PDF + PNG + TIFF)
)

for (s in scripts) {
  message("Running ", s)
  source(file.path(PROJ, "R", s), chdir = TRUE, local = new.env())
}

message("All figures written to ", file.path(PROJ, "figures", "R"))

# -----------------------------------------------------------------------------
# Figure index (script, data source, scientific meaning)
# -----------------------------------------------------------------------------
# 01 class counts              gene_master_table.tsv             locus composition
# 02 sequence features         seq_features_by_class.tsv         promoter syntax axis
# 03 promoter shape            promoter_ext/A_shape_by_class      CAGE peak geometry
# 04 cofactor occupancy        tf_contexts/A_cofactor_occupancy   machinery load
# 05 eta2 taxonomy             promoter_ext/C_class_vs_biotype    taxonomy gain
# 06 TE composition            te_by_class.tsv                    promoter TE content
# 07 TE donor age              te_age_by_class.tsv                donor age
# 08 conservation partial      robustness_partial_corr + phylop    GC-controlled conservation
# 09 phast vs phylop           phylop/phastcons_vs_phylop         cross-track agreement
# 10 substitution              comparative/G_substitution_tests   human-macaque turnover
# 11 TF evidence               motif_enrichment + tf_binding      dual-evidence TFs
# 12 TF pairs                  tf_contexts/A_tf_pairs             motif modules
# 13 TF tissue coupling        tf_contexts/T_tf_class_coupling    tissue TF programmes
# 14 CAGE activity             cage_activity_by_class.tsv         activity/specificity
# 15 ChromHMM                  chromhmm_by_class.tsv              chromatin states
# 16 constraint profiles       constraint_structure/A_*           promoter vs body
# 17 conserved patches         constraint_structure/B_*           patch structure
# 18 TE dual role              constraint_structure/C_*           promoter vs exon
# 19 TT-seq                    ttseq/tests.tsv                    synthesis vs decay
# 20 GTEx composition          tissue_disease/gtex_class_by_tissue tissue composition
# 21 TCGA class DE             tissue_disease/tcga_*              tumour dysregulation
# 22 CGGA                     ccle_cgga/cgga_*                   grade/IDH validation
# 23 CCLE + EPIC               ccle_cgga/ccle_*, methylation      lineage + correlations
# 24 scRNA                    scrna/*                            single-cell validation
# 25 medulloblastoma           cns_cohorts/mb_per_sample          third cohort
# 26 survival/bridge/dev       evo_contexts + cancer_bridge + dev survival + breadth
# 27 age composition           comparative/age_*                  age + method sensitivity
# 28 recent selection          evo_contexts/B_tajima + D_partner  Tajima + partner age
# 29 LTR family                ltr_family/*                       family modules
# 30 mouse taxonomy            mouse_taxonomy/*                   cross-species taxonomy
# 31 cross-species breadth     developmental/cross_species        breadth conservation
# 32 strand controls           strand_controls/*                  parity controls
# 33 strand geometry           strand_parity/*                    length/asymmetry
# 34 ML confusion              ml/A_confusion                     class prediction
# 35 ML ablation               evo_ext/A_ablation + importance    phyloP gain
# 36 perturbation              perturb/*                          DUX4/DNMT validation
# 37 Perturb-seq               perturbseq/*                       TF perturbation
#
# Composites (07_composite_figures.R; rebuilt from the same tables, written as
# PDF, 600-dpi PNG and LZW TIFF at 183 mm width)
# C1 classes + promoters       01-05                              birth-class definition
# C2 TE + conservation + TF    06-13                              repeat-driven origination
# C3 activity + constraint     14-19                              transcription and selection
# C4 tissues + tumours         20-26                              disease context
# C5 evolution + parity        27-33                              age and strand symmetry
# C6 ML + perturbation         34-37                              prediction and causality
