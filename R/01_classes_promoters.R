# =============================================================================
# 01_classes_promoters.R  --  birth classes and promoter architecture
# -----------------------------------------------------------------------------
# Output figures (figures/R/):
#   01_class_counts.pdf/png
#       Data: results/gene_master_table.tsv (module 04).
#       Meaning: number of GENCODE v46 loci in each coordinate-defined class.
#   02_sequence_features.pdf/png
#       Data: results/seq_features_by_class.tsv (module 05).
#       Meaning: GC content, CpG observed/expected, high-CpG fraction,
#                TATA-box fraction by class; the "borrowed vs self-assembled"
#                promoter axis.
#   03_promoter_shape.pdf/png
#       Data: results/promoter_ext/A_shape_by_class.tsv (module 29A,
#             subset = peak_bearing).
#       Meaning: CAGE top-peak width and peak number per promoter; divergent
#                promoters are broad/PCG-like, intergenic/intronic are sharp.
#   04_cofactor_occupancy.pdf/png
#       Data: results/tf_contexts/A_cofactor_occupancy.tsv (module 31A).
#       Meaning: ReMap occupancy of EP300/TFIID/SMARCA4/CTCF across classes.
#   05_eta2_taxonomy.pdf/png
#       Data: results/promoter_ext/C_class_vs_biotype.tsv (module 29C).
#       Meaning: variance explained by the five-class split versus the
#                classic four-way taxonomy (rank-ANOVA eta^2).
# =============================================================================

if (!exists("CLASS_COLORS")) source(file.path(Sys.getenv("LNCPROT_DIR", "."), "R", "00_common.R"))

# ---- 01 class counts --------------------------------------------------------
gm <- read_res("gene_master_table.tsv")
p1 <- gm %>%
  count(cls) %>%
  mutate(cls = fct_cls(cls)) %>%
  ggplot(aes(cls, n, fill = cls)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = comma(n)), vjust = -0.25, size = 2.6) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     labels = comma) +
  labs(x = NULL, y = "Number of loci",
       title = "Birth-class composition of 19,258 GENCODE lncRNAs") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p1, "01_class_counts", 150, 85)

# ---- 02 sequence features ---------------------------------------------------
sf <- read_res("seq_features_by_class.tsv")
sf_long <- sf %>%
  select(cls, gc, cpg_oe, HCP, LCP, TATA, Inr) %>%
  pivot_longer(-cls, names_to = "feature", values_to = "value") %>%
  mutate(cls = fct_cls(cls),
         feature = factor(feature,
                          levels = c("gc", "cpg_oe", "HCP", "LCP", "TATA", "Inr"),
                          labels = c("GC content", "CpG observed/expected",
                                     "High-CpG fraction", "Low-CpG fraction",
                                     "TATA-box fraction", "Initiator fraction")))
p2 <- ggplot(sf_long, aes(cls, value, fill = cls)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.2) +
  facet_wrap(~feature, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Value", title = "Promoter sequence architecture") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p2, "02_sequence_features", 180, 100)

# ---- 03 promoter shape ------------------------------------------------------
sh <- read_res("promoter_ext", "A_shape_by_class.tsv") %>%
  filter(subset == "peak_bearing",
         feature %in% c("top_width", "n_peaks")) %>%
  mutate(cls = fct_cls(cls),
         feature = factor(feature, levels = c("top_width", "n_peaks"),
                          labels = c("Top CAGE-peak width (bp)",
                                     "CAGE peaks per promoter")),
         delta = median_class - median_PCG)
p3 <- ggplot(sh, aes(cls, delta, fill = cls)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  facet_wrap(~feature, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_discrete(labels = CLASS_LABELS) +
  labs(x = NULL, y = "Class median minus PCG median",
       title = "Promoter shape from FANTOM5 CAGE peaks") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p3, "03_promoter_shape", 180, 85)

# ---- 04 cofactor occupancy --------------------------------------------------
cf <- read_res("tf_contexts", "A_cofactor_occupancy.tsv") %>%
  mutate(cls = fct_cls(cls),
         cofactor = factor(cofactor,
                           levels = c("EP300", "CREBBP", "TAF1", "SMARCA4",
                                      "ARID1A", "MED1", "BRD4", "YY1", "CTCF")))
p4 <- ggplot(cf, aes(cofactor, cls, fill = frac_bound)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", frac_bound)), size = 2.1,
            colour = ifelse(cf$frac_bound > 0.7, "white", "black")) +
  scale_fill_gradientn(colours = c("#F7FBFF", "#6BAED6", "#08306B"),
                       limits = c(0, 1), name = "Fraction bound") +
  scale_y_discrete(limits = rev(levels(cf$cls))) +
  labs(x = NULL, y = NULL,
       title = "Cofactor occupancy at promoters (ReMap 2022)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p4, "04_cofactor_occupancy", 170, 95)

# ---- 05 eta^2 taxonomy comparison ------------------------------------------
et <- read_res("promoter_ext", "C_class_vs_biotype.tsv") %>%
  select(feature, eta2_biotype4, eta2_birthclass5) %>%
  pivot_longer(-feature, names_to = "taxonomy", values_to = "eta2") %>%
  mutate(taxonomy = factor(taxonomy,
                           levels = c("eta2_biotype4", "eta2_birthclass5"),
                           labels = c("Classic four-way", "Five birth classes")))
p5 <- ggplot(et, aes(feature, eta2, fill = taxonomy)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("Classic four-way" = "#BBBBBB",
                               "Five birth classes" = "#D55E00"),
                    name = NULL) +
  labs(x = NULL, y = expression(paste("Rank-ANOVA ", eta^2)),
       title = "Variance explained by taxonomy") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top")
save_panel(p5, "05_eta2_taxonomy", 160, 90)
