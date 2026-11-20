# =============================================================================
# 06_ml_perturbation.R  --  machine learning and perturbation validation
# -----------------------------------------------------------------------------
# Output figures:
#   34_ml_confusion.pdf/png
#       Data: results/ml/A_confusion.tsv (module 24).
#       Meaning: out-of-fold confusion matrix for birth-class prediction from
#                promoter-intrinsic features.
#   35_ml_ablation.pdf/png
#       Data: results/evo_ext/A_ablation.tsv,
#             results/evo_ext/A_feature_importance_F2.tsv (module 28A).
#       Meaning: macro-F1 with and without phyloP features, plus seeded
#                permutation importance.
#   36_perturbation.pdf/png
#       Data: results/perturb/perturb_tests.tsv,
#             results/perturb/response_medians.tsv (module 30).
#       Meaning: LTR-derived lncRNA response to DUX4 activation and DNMT
#                inhibition across four public datasets.
#   37_perturbseq.pdf/png
#       Data: results/perturbseq/headline.tsv,
#             class_tests_within_target.tsv (module 39).
#       Meaning: K562 Perturb-seq class-level response to CTCF/MAZ/E2F6/EP300.
# =============================================================================

if (!exists("CLASS_COLORS")) source(file.path(Sys.getenv("LNCPROT_DIR", "."), "R", "00_common.R"))

# ---- 34 ML confusion matrix -------------------------------------------------
cm <- suppressMessages(read_res("ml", "A_confusion.tsv"))
names(cm)[1] <- "true_class"   # first column is unnamed in the TSV
cm <- cm %>%
  pivot_longer(-true_class, names_to = "predicted", values_to = "n") %>%
  mutate(true_class = fct_lnc5(true_class),
         predicted = fct_lnc5(predicted))
cm <- cm %>%
  group_by(true_class) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup()
p34 <- ggplot(cm, aes(predicted, true_class, fill = fraction)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", fraction)), size = 2.2,
            colour = ifelse(cm$fraction > 0.5, "white", "black")) +
  scale_fill_gradient(low = "white", high = "#08306B", name = "Fraction") +
  scale_y_discrete(limits = rev) +
  labs(x = "Predicted class", y = "True class",
       title = "Birth-class prediction from promoter features") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p34, "34_ml_confusion", 120, 100)

# ---- 35 ML ablation --------------------------------------------------------
ab <- read_res("evo_ext", "A_ablation.tsv") %>%
  mutate(feature_set = factor(feature_set, levels = feature_set))
p35a <- ggplot(ab, aes(feature_set, macro_f1, fill = delta_vs_F0 > 0)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.3f", macro_f1)), vjust = -0.4,
            size = 2.4) +
  geom_hline(yintercept = ab$macro_f1[ab$feature_set == "F0_base19"],
             linewidth = 0.35, linetype = 2) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#BBBBBB"),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Macro-F1 (5-fold CV)",
       title = "PhyloP feature ablation") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
im <- read_res("evo_ext", "A_feature_importance_F2.tsv") %>%
  slice_max(abs(delta_f1), n = 15) %>%
  mutate(feature = fct_reorder(feature, delta_f1))
p35b <- ggplot(im, aes(feature, delta_f1, fill = delta_f1 > 0)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#BBBBBB"),
                    guide = "none") +
  labs(x = NULL, y = "Delta macro-F1 after permutation",
       title = "Seeded permutation importance") +
  theme_pub()
p35 <- p35a + p35b + plot_layout(widths = c(1, 1.1))
save_panel(p35, "35_ml_ablation", 195, 90)

# ---- 36 perturbation validation --------------------------------------------
pt <- read_res("perturb", "perturb_tests.tsv") %>%
  filter(contrast == "LTR_lnc_vs_other") %>%
  mutate(dataset_orig = dataset,
         dataset_label = factor(dataset,
                          levels = c("GSE85461_iDUX4_fastq",
                                     "GSE85461_lentiDUX4_14h",
                                     "GSE56787_FSHD_vs_ctrl",
                                     "GSE337536_decitabine",
                                     "GSE305426_DNMTi_vs_DMSO"),
                          labels = c("iDUX4 (GSE85461)",
                                     "lenti-DUX4 14 h",
                                     "FSHD endogenous",
                                     "Decitabine (GSE337536)",
                                     "DAC/AZA (GSE305426)")),
         direction = case_when(
           str_detect(dataset_orig, "GSE85461") ~ "DUX4",
           str_detect(dataset_orig, "FSHD") ~ "Endogenous",
           TRUE ~ "DNMT inhibitor"))
p36 <- ggplot(pt, aes(dataset_label, cliffs, fill = direction)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_text(aes(label = star_lab(pval), vjust = ifelse(cliffs > 0, -0.4, 1.1)),
            size = 2.6) +
  scale_fill_manual(values = c(DUX4 = "#D55E00", Endogenous = "#999999",
                               `DNMT inhibitor` = "#0072B2"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) +
  labs(x = NULL, y = "Cliff's delta (LTR vs other lncRNA)",
       title = "Perturbation validation of LTR-driven induction") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top")
save_panel(p36, "36_perturbation", 180, 90)

# ---- 37 Perturb-seq --------------------------------------------------------
ps <- read_res("perturbseq", "class_tests_within_target.tsv") %>%
  filter(contrast == "divergent_vs_other_lnc") %>%
  mutate(target = fct_reorder(target, cliffs))
p37a <- ggplot(ps, aes(target, cliffs, fill = qval < 0.05)) +
    geom_col(width = 0.6, colour = "black", linewidth = 0.2) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#BBBBBB"),
                      guide = "none") +
    labs(x = NULL, y = "Cliff's delta (divergent lncRNAs)",
         title = "Perturb-seq: divergent response") +
    theme_pub()
hl <- read_res("perturbseq", "headline.tsv") %>%
  mutate(target = factor(target, levels = c("CTCF", "MAZ", "E2F6", "EP300")),
         sig = div_vs_otherlnc_q)
p37b <- ggplot(hl, aes(target, div_vs_otherlnc_cliffs, fill = sig < 0.05)) +
  geom_col(width = 0.6, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#BBBBBB"),
                    guide = "none") +
  labs(x = NULL, y = "Divergent minus other lncRNA",
       title = "Target-level summary") +
  theme_pub()
p37 <- p37a + p37b + plot_layout(widths = c(1, 1))
save_panel(p37, "37_perturbseq", 175, 85)
