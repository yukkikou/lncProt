# =============================================================================
# 05_evolution_parity.R  --  age, selection, LTR families, mouse, strand parity
# -----------------------------------------------------------------------------
# Output figures:
#   27_age_composition.pdf/png
#       Data: results/comparative/E_class_x_age.tsv,
#             age_method_comparison.tsv, age_method_strand_diagnostics.tsv.
#       Meaning: age composition by class, method sensitivity, and the
#                plus/minus detection diagnostics.
#   28_recent_selection.pdf/png
#       Data: results/evo_contexts/B_tajima_tests.tsv,
#             results/comparative/D_partner_age_proxy.tsv.
#       Meaning: Tajima's D by class and partner-PCG conservation.
#   29_ltr_family.pdf/png
#       Data: results/ltr_family/family_by_class.tsv, family_tf_pairs.tsv,
#             family_age.tsv.
#       Meaning: LTR family enrichment by class, family-TF modules and age.
#   30_mouse_taxonomy.pdf/png
#       Data: results/mouse_taxonomy/cross_species_class_comparison.tsv,
#             cross_species_activity_comparison.tsv.
#       Meaning: mouse class composition and activity hierarchy versus human.
#   31_cross_species_breadth.pdf/png
#       Data: results/developmental/cross_species.tsv.
#       Meaning: conservation of developmental expression breadth for
#                syntenic orthologs.
#   32_strand_controls.pdf/png
#       Data: results/strand_controls/S_sources.tsv, T_te_orientation.tsv,
#             cage_version_plus.tsv.
#       Meaning: source balance, TE orientation, and CAGE plus-bias controls.
#   33_strand_geometry.pdf/png
#       Data: results/strand_parity/C_length_parity.tsv,
#             results/strand_parity/E_ltr_asymmetry.tsv,
#             results/strand_parity/E_divergent_tpm.tsv.
#       Meaning: length parity/power-law, LTR sense/antisense symmetry and
#                divergent-pair CAGE correlation.
# =============================================================================

if (!exists("CLASS_COLORS")) source(file.path(Sys.getenv("LNCPROT_DIR", "."), "R", "00_common.R"))

# ---- 27 age composition and method sensitivity ------------------------------
age <- read_res("comparative", "E_class_x_age.tsv") %>%
  pivot_longer(-cls, names_to = "age", values_to = "n") %>%
  mutate(cls = fct_lnc5(cls),
         age = factor(age, levels = rev(AGE_LEVELS),
                      labels = rev(AGE_LABELS)))
p27a <- ggplot(age, aes(cls, n, fill = age)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.15) +
  scale_fill_manual(values = age_pal(), name = NULL) +
  scale_x_discrete(labels = CLASS_LABELS[LNC5]) +
  labs(x = NULL, y = "Loci", title = "Transcript-age composition") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom")
am <- read_res("comparative", "age_method_comparison.tsv") %>%
  filter(cls %in% LNC5, !is.na(human_specific)) %>%
  mutate(cls = fct_lnc5(cls),
         method = factor(method,
                         levels = c("same_strand", "either_strand", "syntenic"),
                         labels = c("Same-strand TSS", "Either-strand TSS",
                                    "Syntenic body")))
p27b <- ggplot(am, aes(method, human_specific, colour = cls, group = cls)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = CLASS_COLORS[LNC5],
                      labels = CLASS_LABELS[LNC5], name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Human-specific fraction",
       title = "Age-assignment method sensitivity") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "bottom")
sd <- read_res("comparative", "age_method_strand_diagnostics.tsv") %>%
  pivot_longer(c(plus, minus), names_to = "strand", values_to = "fraction") %>%
  mutate(method = factor(method,
                         levels = c("same_strand", "either_strand", "syntenic"),
                         labels = c("Same-strand TSS", "Either-strand TSS",
                                    "Syntenic body")),
         species = factor(species, levels = c("macaque", "mouse", "rat")))
p27c <- ggplot(sd, aes(species, fraction, fill = strand)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           colour = "black", linewidth = 0.2) +
  facet_wrap(~method, nrow = 1) +
  scale_fill_manual(values = c(plus = "#0072B2", minus = "#D55E00"),
                    name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Ortholog detection rate",
       title = "Strand diagnostics of age assignment") +
  theme_pub() +
  theme(legend.position = "bottom")
p27 <- (p27a | p27b | p27c) + plot_annotation(tag_levels = "a")
save_panel(p27, "27_age_composition", 220, 95)

# ---- 28 recent selection and partner age -----------------------------------
tj <- read_res("evo_contexts", "B_tajima_tests.tsv") %>%
  filter(feature == "tajD_mean", cls %in% c("PCG", LNC5)) %>%
  mutate(cls = fct_cls(cls), delta = median_class - median_PCG)
p28a <- ggplot(tj, aes(cls, delta, fill = cls)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  labs(x = NULL, y = "Tajima's D minus PCG",
       title = "Recent selection at promoters") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
pt <- read_res("evo_contexts", "D_partner_age_proxy.tsv") %>%
  filter(contrast == "partner_vs_nonpartner_PCG", feature == "phylop_mean")
p28b <- ggplot(pt, aes("Partner PCGs", median_partner)) +
  geom_col(width = 0.5, fill = "#0072B2", colour = "black", linewidth = 0.2) +
  geom_col(aes(x = "Non-partner PCGs", y = median_other), width = 0.5,
           fill = "#BBBBBB", colour = "black", linewidth = 0.2) +
  labs(x = NULL, y = "Median phyloP", title = "Partner-PCG conservation") +
  theme_pub()
p28 <- p28a + p28b + plot_layout(widths = c(1.4, 0.8))
save_panel(p28, "28_recent_selection", 165, 85)

# ---- 29 LTR family modules --------------------------------------------------
lf <- read_res("ltr_family", "family_by_class.tsv") %>%
  filter(cls %in% LNC5) %>%
  mutate(cls = fct_lnc5(cls))
p29a <- ggplot(lf, aes(family, cls, fill = log2OR)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, name = "log2OR") +
  labs(x = NULL, y = NULL, title = "LTR family by class") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5))
ft <- read_res("ltr_family", "family_tf_pairs.tsv") %>%
  group_by(family) %>%
  slice_max(log2OR, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = paste(family, tf, sep = " + "),
         label = fct_reorder(label, log2OR))
p29b <- ggplot(ft, aes(label, log2OR, fill = family)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = rep(c("#D55E00", "#0072B2", "#009E73",
                                   "#CC79A7", "#E69F00", "#56B4E9",
                                   "#999999", "#666666", "#BBBBBB"), 3),
                    guide = "none") +
  labs(x = NULL, y = "Family-TF log2OR",
       title = "Top LTR-family TF modules") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 5.8))
fa <- read_res("ltr_family", "family_age.tsv") %>%
  pivot_longer(-family, names_to = "age", values_to = "n") %>%
  filter(age != "unclassified") %>%
  mutate(age = factor(age, levels = c("ancient", "primate_old",
                                      "primate_young")))
p29c <- ggplot(fa, aes(family, n, fill = age)) +
  geom_col(width = 0.75, colour = "black", linewidth = 0.15) +
  scale_fill_manual(values = c(ancient = "#0072B2",
                               primate_old = "#56B4E9",
                               primate_young = "#D55E00"),
                    name = NULL) +
  coord_flip() +
  labs(x = NULL, y = "TSS-donating elements",
       title = "LTR family age") +
  theme_pub() +
  theme(legend.position = "bottom")
p29 <- p29a | (p29b / p29c)
save_panel(p29, "29_ltr_family", 210, 110)

# ---- 30 mouse taxonomy ------------------------------------------------------
mc <- read_res("mouse_taxonomy", "cross_species_class_comparison.tsv") %>%
  pivot_longer(c(human_fraction, mouse_fraction), names_to = "species",
               values_to = "fraction") %>%
  mutate(cls = fct_cls(cls),
         species = factor(species, levels = c("human_fraction", "mouse_fraction"),
                          labels = c("Human", "Mouse")))
p30a <- ggplot(mc, aes(cls, fraction, fill = species)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c(Human = "#0072B2", Mouse = "#E69F00"),
                    name = NULL) +
  scale_x_discrete(labels = CLASS_LABELS) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Class fraction",
       title = "Class composition in mouse and human") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "top")
ma <- read_res("mouse_taxonomy", "cross_species_activity_comparison.tsv") %>%
  filter(cls %in% c("PCG", LNC5)) %>%
  mutate(cls = fct_cls(cls))
p30b <- ggplot(ma, aes(cls, active_frac, fill = species)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c(human = "#0072B2", mouse = "#E69F00"),
                    labels = c("Human", "Mouse"), name = NULL) +
  scale_x_discrete(labels = CLASS_LABELS) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "CAGE-active fraction",
       title = "Activity hierarchy in mouse and human") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "top")
p30 <- p30a + p30b + plot_layout(widths = c(1, 1.1))
save_panel(p30, "30_mouse_taxonomy", 195, 90)

# ---- 31 cross-species developmental breadth --------------------------------
cs <- read_res("developmental", "cross_species.tsv") %>%
  filter(cls %in% c("PCG", LNC5)) %>%
  mutate(cls = fct_cls(cls),
         species = factor(species, levels = SPECIES_ORDER))
p31 <- ggplot(cs, aes(species, cls, fill = rho_breadth)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n >= 20, sprintf("%.2f", rho_breadth), "")),
            size = 2.1) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, name = "rho") +
  labs(x = NULL, y = NULL,
       title = "Conservation of developmental expression breadth") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p31, "31_cross_species_breadth", 160, 90)

# ---- 32 strand controls -----------------------------------------------------
ss <- read_res("strand_controls", "S_sources.tsv") %>%
  mutate(source = fct_reorder(source, plus_frac))
p32a <- ggplot(ss, aes(source, plus_frac)) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2,
                linewidth = 0.4) +
  geom_hline(yintercept = 0.5, linewidth = 0.35, linetype = 2) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Plus-strand fraction",
       title = "Strand balance by source") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6))
to <- read_res("strand_controls", "T_te_orientation.tsv") %>%
  filter(str_detect(group, "^LTR:|^LTR$")) %>%
  mutate(group = str_remove(group, "^LTR:"),
         group = fct_reorder(group, same_frac))
p32b <- ggplot(to, aes(group, same_frac)) +
  geom_col(width = 0.7, fill = "#D55E00", colour = "black",
           linewidth = 0.2) +
  geom_hline(yintercept = 0.5, linewidth = 0.35, linetype = 2) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Same-strand fraction",
       title = "TE orientation at functional TSS") +
  theme_pub()
cv <- read_res("strand_controls", "cage_version_plus.tsv") %>%
  mutate(set = fct_reorder(set, plus_frac))
p32c <- ggplot(cv, aes(set, plus_frac)) +
  geom_col(width = 0.7, fill = "#0072B2", colour = "black",
           linewidth = 0.2) +
  geom_hline(yintercept = 0.5, linewidth = 0.35, linetype = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 0.6)) +
  labs(x = NULL, y = "Plus-strand fraction",
       title = "CAGE version control") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
p32 <- p32a | (p32b / p32c)
save_panel(p32, "32_strand_controls", 200, 110)

# ---- 33 length parity and transcription asymmetry --------------------------
lp <- read_res("strand_parity", "C_length_parity.tsv") %>%
  filter(feature %in% c("gene", "exon"), cls %in% c("PCG", LNC5)) %>%
  mutate(cls = fct_cls(cls))
p33a <- ggplot(lp, aes(cls, ks_stat, fill = cls)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  facet_wrap(~feature, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  labs(x = NULL, y = "KS statistic (plus vs minus)",
       title = "Length-distribution parity") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
# Most LTR-derived promoters have signal on one strand only, so this panel
# holds tens of values, not thousands. A histogram implied a sample size the
# data does not have; a dot plot with the count stated is honest about it.
la <- read_res("strand_parity", "E_ltr_asymmetry.tsv") %>%
  filter(is.finite(log2_ratio))
p33b <- ggplot(la, aes(x = log2_ratio)) +
  geom_dotplot(binwidth = 0.08, method = "histodot", fill = "#56B4E9",
               colour = "black", stroke = 0.2, dotsize = 0.5) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2) +
  annotate("text", x = Inf, y = Inf, label = sprintf("n = %d", nrow(la)),
           hjust = 1.2, vjust = 1.8, size = 2.2) +
  labs(x = "log2 sense / antisense CAGE", y = "Promoters",
       title = "LTR promoter symmetry") +
  theme_pub()
dt <- read_res("strand_parity", "E_divergent_tpm.tsv") %>%
  filter(!is.na(lnc_tpm), !is.na(partner_tpm))
p33c <- ggplot(dt, aes(lnc_tpm, partner_tpm)) +
  geom_point(size = 0.5, alpha = 0.35, colour = "#E69F00") +
  geom_smooth(method = "lm", se = FALSE, colour = "#0072B2",
              linewidth = 0.7) +
  scale_x_continuous(trans = pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  labs(x = "lncRNA CAGE (mean TPM)", y = "Partner PCG CAGE (mean TPM)",
       title = "Divergent-pair CAGE correlation") +
  theme_pub()
p33 <- p33a | (p33b / p33c)
save_panel(p33, "33_strand_geometry", 210, 110)
