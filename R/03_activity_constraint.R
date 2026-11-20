# =============================================================================
# 03_activity_constraint.R  --  activity, chromatin, promoter-body constraint
# -----------------------------------------------------------------------------
# Output figures:
#   14_cage_activity.pdf/png
#       Data: results/cage_activity_by_class.tsv (module 09).
#       Meaning: CAGE-active fraction, tissue specificity tau, and breadth.
#   15_chromhmm.pdf/png
#       Data: results/chromhmm_by_class.tsv (module 10).
#       Meaning: promoter-like / enhancer-like / bivalent / repressed state
#                fractions at TSS (mean of five Roadmap cell types).
#   16_constraint_profiles.pdf/png
#       Data: results/constraint_structure/A_conservation_tests.tsv (module 43).
#       Meaning: promoter vs gene-body vs exon vs intron phyloP differences
#                relative to PCG.
#   17_conserved_patches.pdf/png
#       Data: results/constraint_structure/B_patch_tests.tsv (module 43).
#       Meaning: conserved-patch fraction and longest patch at promoters.
#   18_te_dual_role.pdf/png
#       Data: results/constraint_structure/C_te_dual_tests.tsv (module 43).
#       Meaning: TE coverage and TE/LTR exonization by class.
#   19_ttseq.pdf/png
#       Data: results/ttseq/tests.tsv (module 35).
#       Meaning: synthesis, decay and synthesis-minus-decay from TT-seq.
# =============================================================================

if (!exists("CLASS_COLORS")) source(file.path(Sys.getenv("LNCPROT_DIR", "."), "R", "00_common.R"))

# ---- 14 CAGE activity -------------------------------------------------------
ac <- read_res("cage_activity_by_class.tsv") %>%
  filter(cls %in% c("PCG", LNC5)) %>%
  mutate(cls = fct_cls(cls))
p14a <- ggplot(ac, aes(cls, active_frac, fill = cls)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "CAGE-active fraction", title = "Promoter activity") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 5.8))
p14b <- ggplot(ac, aes(cls, tau, fill = cls)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  labs(x = NULL, y = "Tissue specificity tau", title = "Specificity") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 5.8))
p14c <- ggplot(ac, aes(cls, breadth, fill = cls)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  labs(x = NULL, y = "CAGE samples > 1 TPM", title = "Expression breadth") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 5.8))
# Do not add `& theme_pub()` here. It re-applies axis.text to every subplot and
# silently resets the per-panel axis.text.x rotation set above.
p14 <- p14a + p14b + p14c +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9))
# Vertical labels need room under each of the three narrow panels.
save_panel(p14, "14_cage_activity", 190, 95)

# ---- 15 ChromHMM states -----------------------------------------------------
ch <- read_res("chromhmm_by_class.tsv") %>%
  filter(cls %in% c("PCG", LNC5)) %>%
  mutate(cls = fct_cls(cls))
ch_mean <- ch %>%
  select(cls, matches("_(prom_like|enh_like|biv_like|repressed)$")) %>%
  pivot_longer(-cls, names_to = "column", values_to = "fraction") %>%
  mutate(state = str_extract(column, "(prom_like|enh_like|biv_like|repressed)$"),
         state = factor(state,
                        levels = c("prom_like", "enh_like", "biv_like",
                                   "repressed"),
                        labels = c("Promoter-like", "Enhancer-like",
                                   "Bivalent", "Repressed"))) %>%
  group_by(cls, state) %>%
  summarise(fraction = mean(fraction, na.rm = TRUE), .groups = "drop")
p15 <- ggplot(ch_mean, aes(cls, fraction, fill = state)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72,
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("Promoter-like" = "#0072B2",
                               "Enhancer-like" = "#E69F00",
                               Bivalent = "#009E73",
                               Repressed = "#999999"),
                    name = NULL) +
  scale_x_discrete(labels = CLASS_LABELS) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Mean fraction of TSS",
       title = "ChromHMM state at TSS (5 cell types)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top")
save_panel(p15, "15_chromhmm", 175, 90)

# ---- 16 promoter vs body vs exon constraint ---------------------------------
cp <- read_res("constraint_structure", "A_conservation_tests.tsv") %>%
  filter(feature %in% c("prom_phylop", "body_phylop", "exon_phylop",
                        "intron_phylop"),
         cls %in% c("PCG", LNC5)) %>%
  mutate(cls = fct_cls(cls),
         feature = factor(feature,
                          levels = c("prom_phylop", "body_phylop",
                                     "exon_phylop", "intron_phylop"),
                          labels = c("Promoter", "Gene body", "Exon", "Intron")),
         delta = median_class - median_PCG)
p16 <- ggplot(cp, aes(cls, delta, fill = cls)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  facet_wrap(~feature, nrow = 1) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  labs(x = NULL, y = "Class median minus PCG median (phyloP)",
       title = "Constraint is concentrated at the promoter") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(p16, "16_constraint_profiles", 190, 85)

# ---- 17 conserved patches ---------------------------------------------------
pt <- read_res("constraint_structure", "B_patch_tests.tsv") %>%
  filter(feature %in% c("prom_patch_frac_phylop", "prom_patch_max_phylop"),
         cls %in% c("PCG", LNC5)) %>%
  mutate(cls = fct_cls(cls),
         feature = factor(feature,
                          levels = c("prom_patch_frac_phylop",
                                     "prom_patch_max_phylop"),
                          labels = c("Fraction of promoter in patches",
                                     "Longest conserved patch (bp)")),
         delta = median_class - median_PCG)
p17 <- ggplot(pt, aes(cls, delta, fill = cls)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  facet_wrap(~feature, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  labs(x = NULL, y = "Class median minus PCG median",
       title = "Conserved-patch structure (phyloP >= 2)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(p17, "17_conserved_patches", 185, 85)

# ---- 18 TE promoter vs exon -------------------------------------------------
td <- read_res("constraint_structure", "C_te_dual_tests.tsv") %>%
  filter(feature %in% c("te_cov", "te_exon_frac", "ltr_exon_frac"),
         cls %in% c("PCG", LNC5) | str_detect(cls, "LTR")) %>%
  mutate(cls = fct_cls(cls),
         feature = factor(feature,
                          levels = c("te_cov", "te_exon_frac",
                                     "ltr_exon_frac"),
                          labels = c("Promoter TE coverage",
                                     "TE-derived exon fraction",
                                     "LTR-derived exon fraction")),
         delta = median_class - median_PCG)
p18 <- ggplot(td, aes(cls, delta, fill = cls)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  facet_wrap(~feature, nrow = 1) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  labs(x = NULL, y = "Class median minus PCG median",
       title = "TE promoter donation versus exonization") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(p18, "18_te_dual_role", 190, 85)

# ---- 19 TT-seq synthesis and decay -----------------------------------------
tt <- read_res("ttseq", "tests.tsv") %>%
  filter(feature %in% c("log_syn", "log_decay", "log_ratio"),
         cls %in% c("PCG", LNC5)) %>%
  mutate(cls = fct_cls(cls),
         feature = factor(feature,
                          levels = c("log_syn", "log_decay", "log_ratio"),
                          labels = c("log10 synthesis", "log10 decay",
                                     "log10 synthesis minus decay")),
         delta = median_class - median_PCG)
p19 <- ggplot(tt, aes(cls, delta, fill = cls)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  facet_wrap(~feature, nrow = 1) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  labs(x = NULL, y = "Class median minus PCG median",
       title = "TT-seq: synthesis versus decay") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(p19, "19_ttseq", 190, 85)
