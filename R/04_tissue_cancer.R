# =============================================================================
# 04_tissue_cancer.R  --  GTEx, TCGA and validation cohorts
# -----------------------------------------------------------------------------
# Output figures:
#   20_gtex_composition.pdf/png
#       Data: results/tissue_disease/tables/gtex_class_by_tissue.tsv (13c).
#       Meaning: birth-class share of expressed lncRNAs across 49 GTEx tissues.
#   21_tcga_class_de.pdf/png
#       Data: results/tissue_disease/tables/tcga_class_de_per_cancer.tsv,
#             tcga_ltr_derepression.tsv, tcga_class_stage.tsv (13d).
#       Meaning: fraction up-regulated per cancer, LTR derepression, stage.
#   22_cgga.pdf/png
#       Data: results/ccle_cgga/tables/cgga_class_by_group.tsv,
#             cgga_ltr_grade_idh_tests.tsv (19).
#       Meaning: class detection across WHO grade/IDH and LTR-relative delta.
#   23_ccle_epic.pdf/png
#       Data: results/ccle_cgga/tables/ccle_class_by_lineage.tsv,
#             ccle_dux4_vs_ltr.tsv, methylation_vs_ltr.tsv (19-20).
#       Meaning: lineage hierarchy and DUX4/methylation correlations.
#   24_scrna.pdf/png
#       Data: results/scrna/class_by_celltype.tsv,
#             divergent_co_detection.tsv (23).
#       Meaning: class detection by cell type and divergent-pair co-detection.
#   25_medulloblastoma.pdf/png
#       Data: results/cns_cohorts/mb_per_sample.tsv (26).
#       Meaning: per-sample class fractions and LTR-relative activation.
#   26_survival_bridge_dev.pdf/png
#       Data: results/evo_contexts/A_survival_meta.tsv,
#             results/cancer_bridge/C_bridge_tests.tsv,
#             results/developmental/human_class_tests.tsv.
#       Meaning: survival meta, cancer-association odds ratios, and human
#                developmental breadth by class.
# =============================================================================

if (!exists("CLASS_COLORS")) source(file.path(Sys.getenv("LNCPROT_DIR", "."), "R", "00_common.R"))

# ---- 20 GTEx composition and specificity -----------------------------------
gt <- read_res("tissue_disease", "tables", "gtex_class_by_tissue.tsv")
# Only the five birth-class shares partition the expressed set; share_LTR is an
# overlapping annotation and must not enter the stack or totals exceed 1.
SHARE_COLS <- paste0("share_", LNC5)
gt_long <- gt %>%
  select(tissue, all_of(SHARE_COLS), n_lnc_expressed) %>%
  pivot_longer(all_of(SHARE_COLS), names_to = "cls", values_to = "share") %>%
  mutate(cls = str_remove(cls, "^share_"),
         cls = factor(cls, levels = LNC5),
         tissue = fct_reorder(tissue, n_lnc_expressed))
p20a <- ggplot(gt_long, aes(tissue, share, fill = cls)) +
  geom_col(width = 0.8, colour = "black", linewidth = 0.1) +
  coord_flip() +
  scale_fill_manual(values = CLASS_COLORS[LNC5],
                    labels = CLASS_LABELS[LNC5], name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Share of expressed lncRNAs",
       title = "GTEx tissue composition") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 5.2),
        legend.position = "bottom")
sp <- read_res("tissue_disease", "tables", "gtex_class_specificity.tsv") %>%
  filter(cls %in% c("PCG", LNC5)) %>%
  mutate(cls = fct_cls(cls))
p20b <- ggplot(sp, aes(cls, tau_median, fill = cls)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  labs(x = NULL, y = "Median GTEx tau",
       title = "Tissue specificity") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
p20 <- p20a + p20b + plot_layout(widths = c(1.5, 1))
save_panel(p20, "20_gtex_composition", 190, 110)

# ---- 21 TCGA class DE / LTR / stage ----------------------------------------
# DE was only run in cohorts with at least five adjacent normals, so frac_up is
# undefined for the remaining cohorts and fct_reorder() cannot rank them. Drop
# those rows explicitly instead of letting NaN propagate into the ordering.
de <- read_res("tissue_disease", "tables", "tcga_class_de_per_cancer.tsv") %>%
  filter(cls %in% LNC5, n_normal >= 5) %>%
  mutate(cls = fct_lnc5(cls),
         cancer = fct_reorder(cancer, frac_up, mean, na.rm = TRUE))
p21a <- ggplot(de, aes(cancer, cls, fill = frac_up)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "#08306B",
                      name = "Fraction up") +
  labs(x = NULL, y = NULL, title = "TCGA class up-regulation") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))
ltr <- read_res("tissue_disease", "tables", "tcga_ltr_derepression.tsv") %>%
  mutate(cancer = fct_reorder(cancer, diff))
p21b <- ggplot(ltr, aes(diff, cancer, fill = diff > 0)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#009E73"),
                    guide = "none") +
  labs(x = "LTR-derived minus other lncRNA (log2TPM)", y = NULL,
       title = "LTR derepression") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 5.5))
stg <- read_res("tissue_disease", "tables", "tcga_class_stage.tsv") %>%
  filter(cls %in% LNC5) %>%
  mutate(cls = fct_lnc5(cls))
p21c <- ggplot(stg, aes(cancer, cls, fill = median_delta_advanced)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, name = "Delta") +
  labs(x = NULL, y = NULL, title = "Advanced vs early stage") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))
p21 <- p21a + p21b + p21c + plot_layout(widths = c(1.3, 0.9, 1.1))
save_panel(p21, "21_tcga_class_de", 200, 100)

# ---- 22 CGGA grade and IDH --------------------------------------------------
cg <- read_res("ccle_cgga", "tables", "cgga_class_by_group.tsv") %>%
  filter(group == "Grade", cls %in% c("LTR_derived_lnc", LNC5)) %>%
  mutate(cls = case_when(cls == "LTR_derived_lnc" ~ "LTR_lnc", TRUE ~ cls),
         level = factor(level, levels = c("WHO II", "WHO III", "WHO IV")))
p22a <- ggplot(cg, aes(level, frac_expr, colour = cls, group = cls)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c(CLASS_COLORS[LNC5],
                                 LTR_lnc = "#D55E00"), name = NULL) +
  labs(x = NULL, y = "Fraction of class genes detected",
       title = "CGGA glioma grade") +
  theme_pub() +
  theme(legend.position = "bottom", legend.text = element_text(size = 6))
ct <- read_res("ccle_cgga", "tables", "cgga_ltr_grade_idh_tests.tsv")
p22b <- ggplot(ct, aes(contrast, median_A - median_B)) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  coord_flip() +
  labs(x = NULL, y = "Median LTR delta difference",
       title = "Grade / IDH contrasts") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6))
p22 <- p22a + p22b + plot_layout(widths = c(1.2, 1))
save_panel(p22, "22_cgga", 180, 95)

# ---- 23 CCLE lineage and EPIC/CCLE correlations ----------------------------
cl <- read_res("ccle_cgga", "tables", "ccle_class_by_lineage.tsv") %>%
  filter(cls %in% LNC5) %>%
  mutate(cls = fct_lnc5(cls))
top_lin <- cl %>% distinct(lineage, n_lines) %>% slice_max(n_lines, n = 15) %>%
  pull(lineage)
cl <- cl %>% filter(lineage %in% top_lin) %>%
  mutate(lineage = fct_reorder(lineage, n_lines))
p23a <- ggplot(cl, aes(lineage, cls, fill = frac_expr)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "#08306B",
                      name = "Fraction") +
  labs(x = NULL, y = NULL, title = "CCLE lineage hierarchy") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))
corr_epic <- read_res("ccle_cgga", "tables", "methylation_vs_ltr.tsv") %>%
  mutate(dataset = "EPIC methylation", estimate = rho) %>%
  select(dataset, estimate)
corr_ccle <- read_res("ccle_cgga", "tables", "ccle_dux4_vs_ltr.tsv") %>%
  mutate(dataset = "CCLE DUX4", estimate = rho) %>%
  select(dataset, estimate)
p23b <- bind_rows(corr_epic, corr_ccle) %>%
  ggplot(aes(estimate, dataset, fill = estimate > 0)) +
  geom_col(width = 0.55, colour = "black", linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"),
                    guide = "none") +
  labs(x = "Spearman rho", y = NULL, title = "LTR correlates") +
  theme_pub()
p23 <- p23a + p23b + plot_layout(widths = c(1.5, 1))
save_panel(p23, "23_ccle_epic", 180, 95)

# ---- 24 scRNA ---------------------------------------------------------------
sc <- read_res("scrna", "class_by_celltype.tsv") %>%
  filter(cls %in% c("PCG", LNC5)) %>%
  mutate(cls = fct_cls(cls))
p24a <- ggplot(sc, aes(cls, cell_type, fill = det_rate)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "#08306B",
                      name = "Detection") +
  labs(x = NULL, y = NULL, title = "scRNA class detection") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5))
cd <- read_res("scrna", "divergent_co_detection.tsv")
p24b <- ggplot(cd, aes(ratio)) +
  geom_histogram(bins = 35, fill = "#56B4E9", colour = "black",
                 linewidth = 0.2) +
  geom_vline(xintercept = 1, linewidth = 0.35, linetype = 2) +
  labs(x = "Observed / expected co-detection", y = "Pairs",
       title = "Divergent-pair co-detection") +
  theme_pub()
p24 <- p24a + p24b + plot_layout(widths = c(1.2, 1))
save_panel(p24, "24_scrna", 175, 95)

# ---- 25 medulloblastoma -----------------------------------------------------
mb <- read_res("cns_cohorts", "mb_per_sample.tsv") %>%
  select(sample, matches("_frac$"), ltr_rel) %>%
  pivot_longer(-c(sample, ltr_rel), names_to = "cls", values_to = "fraction") %>%
  mutate(cls = str_remove(cls, "_frac$")) %>%
  filter(cls %in% LNC5) %>%
  mutate(cls = factor(cls, levels = LNC5))
p25a <- ggplot(mb, aes(cls, fraction, fill = cls)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
  scale_fill_manual(values = CLASS_COLORS[LNC5], guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS[LNC5]) +
  labs(x = NULL, y = "Per-sample class fraction",
       title = "Medulloblastoma cohort") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
p25b <- ggplot(read_res("cns_cohorts", "mb_per_sample.tsv"),
               aes(ltr_rel)) +
  geom_histogram(bins = 30, fill = "#D55E00", colour = "black",
                 linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2) +
  labs(x = "LTR-relative activation", y = "Samples",
       title = "LTR activation") +
  theme_pub()
p25 <- p25a + p25b + plot_layout(widths = c(1.2, 1))
save_panel(p25, "25_medulloblastoma", 175, 90)

# ---- 26 survival, cancer bridge, developmental -----------------------------
sv <- read_res("evo_contexts", "A_survival_meta.tsv") %>%
  mutate(feature = fct_cls(feature))
p26a <- ggplot(sv, aes(feature, z_meta, fill = z_meta > 0)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"),
                    guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  labs(x = NULL, y = "Stouffer z", title = "TCGA survival meta") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
br <- read_res("cancer_bridge", "C_bridge_tests.tsv") %>%
  filter(str_detect(contrast, "^OR_")) %>%
  mutate(cls = str_remove(contrast, "^OR_"),
         cls = fct_cls(cls))
p26b <- ggplot(br, aes(cls, stat, fill = cls)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 1, linewidth = 0.35, linetype = 2) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  labs(x = NULL, y = "Odds ratio (PC-like vs rest)",
       title = "Cancer-association bridge") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
dv <- read_res("developmental", "human_class_tests.tsv") %>%
  filter(feature %in% c("breadth_frac", "tau"), cls %in% c("PCG", LNC5)) %>%
  mutate(cls = fct_cls(cls),
         feature = factor(feature, levels = c("breadth_frac", "tau"),
                          labels = c("Developmental breadth", "Tau")))
p26c <- ggplot(dv, aes(cls, median_class, fill = cls)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  facet_wrap(~feature, scales = "free_y") +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  labs(x = NULL, y = "Median", title = "Human developmental series") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
p26 <- (p26a | p26b | p26c) + plot_annotation(tag_levels = "a")
save_panel(p26, "26_survival_bridge_dev", 200, 85)
