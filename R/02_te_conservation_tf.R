# =============================================================================
# 02_te_conservation_tf.R  --  TE content, conservation, TF evidence
# -----------------------------------------------------------------------------
# Output figures:
#   06_te_composition.pdf/png
#       Data: results/te_by_class.tsv (module 06).
#       Meaning: promoter coverage by LTR/LINE/SINE/DNA and any TE.
#   07_te_age.pdf/png
#       Data: results/te_age_by_class.tsv (module 12b).
#       Meaning: age composition of TSS-donating TEs (young vs ancient).
#   08_conservation_partial.pdf/png
#       Data: results/robustness_partial_corr.tsv and
#             results/phylop/partial_corr_phylop.tsv (modules 12d/27).
#       Meaning: class-conservation association before and after GC control,
#                for phastCons470 and phyloP241.
#   09_phast_vs_phylop.pdf/png
#       Data: results/phylop/phastcons_vs_phylop.tsv (module 27).
#       Meaning: class-level phastCons versus phyloP, showing cross-track
#                agreement.
#   10_substitution.pdf/png
#       Data: results/comparative/G_substitution_tests.tsv (module 34G).
#       Meaning: human-macaque substitution excess by class, all sites and
#                non-CpG sites.
#   11_tf_evidence.pdf/png
#       Data: results/motif_enrichment.tsv, results/tf_binding_enrichment.tsv
#             (modules 07-08).
#       Meaning: top motif and ChIP enrichment per class (dual evidence).
#   12_tf_pairs.pdf/png
#       Data: results/tf_contexts/A_tf_pairs.tsv (module 31A).
#       Meaning: strongest pairwise motif co-occurrence modules per class.
#   13_tf_tissue_coupling.pdf/png
#       Data: results/tf_contexts/T_tf_class_coupling.tsv (module 31T).
#       Meaning: TF expression versus class activity across GTEx tissues.
# =============================================================================

if (!exists("CLASS_COLORS")) source(file.path(Sys.getenv("LNCPROT_DIR", "."), "R", "00_common.R"))

# ---- 06 TE composition ------------------------------------------------------
te <- read_res("te_by_class.tsv") %>%
  select(cls, LTR, LINE, SINE, DNA) %>%
  pivot_longer(-cls, names_to = "order", values_to = "fraction") %>%
  mutate(cls = fct_cls(cls),
         order = factor(order, levels = TE_LEVELS))
p6 <- ggplot(te, aes(cls, fraction, fill = order)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.15) +
  scale_fill_manual(values = TE_COLORS, name = "TE order") +
  scale_x_discrete(labels = CLASS_LABELS) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Fraction of promoter covered",
       title = "TE composition of promoter windows") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p6, "06_te_composition", 165, 90)

# ---- 07 TE donor age --------------------------------------------------------
ta <- read_res("te_age_by_class.tsv") %>%
  select(cls, frac_young, frac_ancient) %>%
  pivot_longer(-cls, names_to = "age", values_to = "fraction") %>%
  mutate(cls = fct_cls(cls),
         age = factor(age, levels = c("frac_ancient", "frac_young"),
                      labels = c("Ancient", "Primate-young")))
p7 <- ggplot(ta, aes(cls, fraction, fill = age)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c(Ancient = "#0072B2",
                               "Primate-young" = "#D55E00"),
                    name = NULL) +
  scale_x_discrete(labels = CLASS_LABELS) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Fraction of TSS-donating TEs",
       title = "Age of TE donors") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top")
save_panel(p7, "07_te_age", 165, 85)

# ---- 08 conservation before/after GC control --------------------------------
pa <- read_res("robustness_partial_corr.tsv") %>%
  mutate(track = "phastCons 470-way")
pl <- read_res("phylop", "partial_corr_phylop.tsv") %>%
  mutate(track = "phyloP 241-mammal")
pc <- bind_rows(pa, pl) %>%
  filter(cls %in% LNC5) %>%
  mutate(cls = fct_lnc5(cls)) %>%
  select(cls, track, rho_raw, rho_GC_controlled) %>%
  pivot_longer(c(rho_raw, rho_GC_controlled),
               names_to = "control", values_to = "rho") %>%
  mutate(control = factor(control,
                          levels = c("rho_raw", "rho_GC_controlled"),
                          labels = c("Raw", "GC-controlled")))
p8 <- ggplot(pc, aes(cls, rho, fill = control)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  facet_wrap(~track, nrow = 1) +
  scale_fill_manual(values = c(Raw = "#0072B2",
                               "GC-controlled" = "#56B4E9"),
                    name = NULL) +
  scale_x_discrete(labels = CLASS_LABELS[LNC5]) +
  labs(x = NULL, y = "Spearman rho (class vs conservation)",
       title = "Conservation gradients before and after GC control") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top")
save_panel(p8, "08_conservation_partial", 175, 85)

# ---- 09 phastCons vs phyloP -------------------------------------------------
pp <- read_res("phylop", "phastcons_vs_phylop.tsv") %>%
  filter(cls %in% c("PCG", LNC5)) %>%
  mutate(cls = fct_cls(cls))
p9 <- ggplot(pp, aes(phast_mean, phylop_mean, colour = cls)) +
  geom_point(size = 2.2) +
  geom_text(aes(label = CLASS_LABELS[cls]), nudge_y = 0.035,
            size = 2.2, show.legend = FALSE) +
  scale_colour_manual(values = CLASS_COLORS, name = NULL) +
  labs(x = "Mean phastCons (470-way)", y = "Mean phyloP (241-mammal)",
       title = "Cross-track conservation of promoter classes") +
  theme_pub() +
  theme(legend.position = "none")
save_panel(p9, "09_phast_vs_phylop", 110, 90)

# ---- 10 substitution excess -------------------------------------------------
sb <- read_res("comparative", "G_substitution_tests.tsv") %>%
  filter(feature %in% c("rate_all", "rate_non"),
         cls %in% c("PCG", LNC5) | str_detect(cls, "LTR")) %>%
  mutate(cls = fct_cls(cls),
         feature = factor(feature, levels = c("rate_all", "rate_non"),
                          labels = c("All sites", "Non-CpG sites")),
         delta = median_class - median_PCG)
p10 <- ggplot(sb, aes(cls, delta, fill = cls)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  facet_wrap(~feature, ncol = 1) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  labs(x = NULL, y = "Class median minus PCG median (substitutions/site)",
       title = "Human-macaque promoter turnover") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p10, "10_substitution", 165, 105)

# ---- 11 TF motif and binding evidence --------------------------------------
# Undefined log2OR (zero cells) cannot enter a ranking scale, so drop them here.
me <- read_res("motif_enrichment.tsv") %>%
  filter(cls %in% LNC5, !is.na(log2OR)) %>%
  mutate(cls = fct_lnc5(cls))
be <- read_res("tf_binding_enrichment.tsv") %>%
  filter(cls %in% LNC5, !is.na(log2OR)) %>%
  mutate(cls = fct_lnc5(cls))
top_tfs <- bind_rows(me, be) %>%
  group_by(cls) %>%
  slice_max(log2OR, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  pull(tf) %>%
  unique()
ev <- bind_rows(
  me %>% filter(tf %in% top_tfs) %>% mutate(assay = "Motif"),
  be %>% filter(tf %in% top_tfs) %>% mutate(assay = "Binding")
) %>%
  mutate(assay = factor(assay, levels = c("Motif", "Binding")),
         tf = fct_reorder(tf, log2OR, .fun = max))
p11 <- ggplot(ev, aes(tf, cls, fill = log2OR)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  facet_wrap(~assay, nrow = 1) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0, name = "log2OR") +
  scale_y_discrete(limits = rev(levels(ev$cls))) +
  labs(x = NULL, y = NULL, title = "Top dual-evidence TFs by class") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.5))
save_panel(p11, "11_tf_evidence", 190, 85)

# ---- 12 TF pair modules -----------------------------------------------------
tp <- read_res("tf_contexts", "A_tf_pairs.tsv") %>%
  filter(cls %in% LNC5) %>%
  mutate(cls = fct_lnc5(cls)) %>%
  group_by(cls) %>%
  slice_max(orr, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(pair = paste(tf_a, tf_b, sep = " + "),
         pair = fct_reorder(pair, orr))
p12 <- ggplot(tp, aes(pair, orr, fill = cls)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  facet_wrap(~cls, scales = "free_y", nrow = 1) +
  coord_flip() +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  labs(x = NULL, y = "Motif co-occurrence log2OR",
       title = "Top TF-pair modules per class") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 5.8))
save_panel(p12, "12_tf_pairs", 190, 100)

# ---- 13 TF-class tissue coupling -------------------------------------------
tc <- read_res("tf_contexts", "T_tf_class_coupling.tsv") %>%
  filter(scope == "gtex_tissue", cls %in% LNC5, fdr < 0.05) %>%
  mutate(cls = fct_lnc5(cls)) %>%
  group_by(cls) %>%
  slice_max(abs(rho), n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = paste(tf, cls, sep = " / "),
         label = fct_reorder(label, rho))
p13 <- ggplot(tc, aes(label, rho, fill = cls)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  labs(x = NULL, y = "Spearman rho (TF expression vs class activity)",
       title = "Strongest TF-class couplings across GTEx tissues") +
  theme_pub()
save_panel(p13, "13_tf_tissue_coupling", 150, 100)
