# =============================================================================
# 07_composite_figures.R  --  merge the 37 single panels into 6 composite figures
# -----------------------------------------------------------------------------
# Produces six multi-panel figures (PDF, PNG 600 dpi, TIFF 600 dpi LZW):
#   C1  birth classes and promoter architecture      (panels 01-05)
#   C2  TE content, conservation and TF evidence     (panels 06-13)
#   C3  activity, chromatin and constraint           (panels 14-19)
#   C4  tissues, tumours and validation cohorts      (panels 20-26)
#   C5  evolution and strand parity                  (panels 27-33)
#   C6  machine learning and perturbation            (panels 34-37)
#
# Each composite is assembled from the ggplot objects built by the individual
# scripts, so every panel keeps its own data source comment there.
# Run after 00_common.R; requires patchwork.
# =============================================================================

if (!exists("CLASS_COLORS")) source(file.path(Sys.getenv("LNCPROT_DIR", "."), "R", "00_common.R"))

save_composite <- function(plot, name, width = 183, height = 240) {
  w <- width / 25.4
  h <- height / 25.4
  ggsave(file.path(FIGDIR, paste0(name, ".pdf")), plot,
         width = w, height = h, units = "in")
  ggsave(file.path(FIGDIR, paste0(name, ".png")), plot,
         width = w, height = h, units = "in", dpi = 600)
  # TIFF with LZW compression, suitable for journal submission
  ggsave(file.path(FIGDIR, paste0(name, ".tif")), plot,
         width = w, height = h, units = "in", dpi = 600,
         compression = "lzw", device = "tiff")
  invisible(plot)
}

# -----------------------------------------------------------------------------
# Helper: rebuild a panel by evaluating a fragment of a script in a sandbox.
# Simpler and more robust than returning objects across files: each builder
# below re-creates the ggplot from the same tables with the same conventions.
# -----------------------------------------------------------------------------

# ---- C1: classes and promoter architecture ---------------------------------
panel_c1 <- function() {
  gm <- read_res("gene_master_table.tsv")
  p1 <- gm %>% count(cls) %>% mutate(cls = fct_cls(cls)) %>%
    ggplot(aes(cls, n, fill = cls)) +
    geom_col(width = 0.72, colour = "black", linewidth = 0.2) +
    scale_fill_manual(values = CLASS_COLORS, guide = "none") +
    scale_x_discrete(labels = CLASS_LABELS) +
    scale_y_continuous(labels = comma) +
    labs(x = NULL, y = "Loci", title = "Class composition") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))

  sf <- read_res("seq_features_by_class.tsv") %>%
    select(cls, gc, cpg_oe, HCP, TATA) %>%
    pivot_longer(-cls, names_to = "feature", values_to = "value") %>%
    mutate(cls = fct_cls(cls),
           feature = factor(feature, levels = c("gc", "cpg_oe", "HCP", "TATA"),
                            labels = c("GC", "CpG o/e", "High-CpG", "TATA")))
  p2 <- ggplot(sf, aes(cls, value, fill = cls)) +
    geom_col(width = 0.72, colour = "black", linewidth = 0.15) +
    facet_wrap(~feature, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = CLASS_COLORS, guide = "none") +
    scale_x_discrete(labels = CLASS_LABELS) +
    labs(x = NULL, y = "Value", title = "Promoter sequence") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4),
                        axis.text.x = element_text(angle = 90, hjust = 1,
                                                   vjust = 0.5, size = 3.4))

  sh <- read_res("promoter_ext", "A_shape_by_class.tsv") %>%
    filter(subset == "peak_bearing", feature == "top_width") %>%
    mutate(cls = fct_cls(cls))
  p3 <- ggplot(sh, aes(cls, median_class, fill = cls)) +
    geom_col(width = 0.72, colour = "black", linewidth = 0.15) +
    geom_hline(yintercept = sh$median_PCG[1], linetype = 2, linewidth = 0.35) +
    scale_fill_manual(values = CLASS_COLORS, guide = "none") +
    scale_x_discrete(labels = CLASS_LABELS) +
    labs(x = NULL, y = "Top peak width (bp)", title = "Promoter shape") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))

  # Keep five representative cofactors. Dropping rows first avoids factor()
  # turning the other four into an unlabelled NA column at the panel edge.
  COF5 <- c("EP300", "TAF1", "SMARCA4", "CTCF", "YY1")
  cf <- read_res("tf_contexts", "A_cofactor_occupancy.tsv") %>%
    filter(cofactor %in% COF5) %>%
    mutate(cls = fct_cls(cls), cofactor = factor(cofactor, levels = COF5))
  p4 <- ggplot(cf, aes(cofactor, cls, fill = frac_bound)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", frac_bound)), size = 2.2,
              colour = if_else(cf$frac_bound > 0.7, "white", "black")) +
    scale_fill_gradientn(colours = c("#F7FBFF", "#6BAED6", "#08306B"),
                         limits = c(0, 1), name = "Bound") +
    scale_y_discrete(limits = rev(levels(droplevels(cf$cls))),
                     labels = CLASS_LABELS, drop = TRUE) +
    labs(x = NULL, y = NULL, title = "Cofactor occupancy") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))

  et <- read_res("promoter_ext", "C_class_vs_biotype.tsv") %>%
    select(feature, eta2_biotype4, eta2_birthclass5) %>%
    pivot_longer(-feature, names_to = "taxonomy", values_to = "eta2") %>%
    mutate(taxonomy = factor(taxonomy,
             levels = c("eta2_biotype4", "eta2_birthclass5"),
             labels = c("Four-way", "Birth classes")))
  p5 <- ggplot(et, aes(feature, eta2, fill = taxonomy)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7,
             colour = "black", linewidth = 0.15) +
    scale_fill_manual(values = c(`Four-way` = "#BBBBBB",
                                 `Birth classes` = "#D55E00"), name = NULL) +
    labs(x = NULL, y = expression(eta^2), title = "Taxonomy gain") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2),
                        legend.position = "top")

  ((p1 | p2) / (p3 | p4 | p5)) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(tag_levels = "a", title = "C1  Birth classes and promoters") + theme(plot.margin = margin(4, 6, 12, 4))
}

# ---- C2: TE, conservation, TF evidence -------------------------------------
panel_c2 <- function() {
  te <- read_res("te_by_class.tsv") %>%
    select(cls, LTR, LINE, SINE, DNA) %>%
    pivot_longer(-cls, names_to = "order", values_to = "fraction") %>%
    mutate(cls = fct_cls(cls), order = factor(order, levels = TE_LEVELS))
  p6 <- ggplot(te, aes(cls, fraction, fill = order)) +
    geom_col(width = 0.72, colour = "black", linewidth = 0.15) +
    scale_fill_manual(values = TE_COLORS, name = "TE") +
    scale_x_discrete(labels = CLASS_LABELS) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Fraction", title = "TE composition") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2),
                        legend.position = "top")

  ta <- read_res("te_age_by_class.tsv") %>%
    select(cls, frac_young, frac_ancient) %>%
    pivot_longer(-cls, names_to = "age", values_to = "fraction") %>%
    mutate(cls = fct_cls(cls),
           age = factor(age, levels = c("frac_ancient", "frac_young"),
                        labels = c("Ancient", "Young")))
  p7 <- ggplot(ta, aes(cls, fraction, fill = age)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7,
             colour = "black", linewidth = 0.15) +
    scale_fill_manual(values = c(Ancient = "#0072B2", Young = "#D55E00"),
                      name = NULL) +
    scale_x_discrete(labels = CLASS_LABELS) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Fraction", title = "TE donor age") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2),
                        legend.position = "top")

  pc <- bind_rows(
    read_res("robustness_partial_corr.tsv") %>% mutate(track = "phastCons"),
    read_res("phylop", "partial_corr_phylop.tsv") %>% mutate(track = "phyloP")
  ) %>% filter(cls %in% LNC5) %>%
    pivot_longer(c(rho_raw, rho_GC_controlled),
                 names_to = "control", values_to = "rho") %>%
    mutate(cls = fct_lnc5(cls),
           control = factor(control, levels = c("rho_raw", "rho_GC_controlled"),
                            labels = c("Raw", "GC-ctrl")))
  p8 <- ggplot(pc, aes(cls, rho, fill = control)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7,
             colour = "black", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    facet_wrap(~track, nrow = 1) +
    scale_fill_manual(values = c(Raw = "#0072B2", `GC-ctrl` = "#56B4E9"),
                      name = NULL) +
    scale_x_discrete(labels = CLASS_LABELS[LNC5]) +
    labs(x = NULL, y = "Spearman rho", title = "Conservation vs class") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2),
                        legend.position = "top")

  sb <- read_res("comparative", "G_substitution_tests.tsv") %>%
    filter(feature %in% c("rate_all", "rate_non"), cls %in% c("PCG", LNC5)) %>%
    mutate(cls = fct_cls(cls),
           feature = factor(feature, levels = c("rate_all", "rate_non"),
                            labels = c("All", "Non-CpG")),
           delta = median_class - median_PCG)
  p10 <- ggplot(sb, aes(cls, delta, fill = cls)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    facet_wrap(~feature, nrow = 1) +
    scale_fill_manual(values = CLASS_COLORS, guide = "none") +
    scale_x_discrete(labels = CLASS_LABELS) +
    labs(x = NULL, y = "Subst./site delta", title = "Human-macaque turnover") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))

  # Motif and binding log2OR live on very different ranges, so the two assays
  # get their own colour scale; a shared scale washes the motif panel white.
  me <- read_res("motif_enrichment.tsv") %>% filter(cls %in% LNC5, !is.na(log2OR))
  be <- read_res("tf_binding_enrichment.tsv") %>% filter(cls %in% LNC5, !is.na(log2OR))
  pick_top <- function(d, n = 10) {
    d %>% group_by(tf) %>% summarise(m = max(log2OR), .groups = "drop") %>%
      slice_max(m, n = n, with_ties = FALSE) %>% pull(tf)
  }
  tf_heat <- function(d, assay) {
    d <- d %>% filter(tf %in% pick_top(d)) %>%
      mutate(cls = fct_lnc5(cls),
             tf = fct_reorder(tf, log2OR, .fun = max, .na_rm = TRUE))
    ggplot(d, aes(tf, cls, fill = log2OR)) +
      geom_tile(colour = "white", linewidth = 0.25) +
      scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                           midpoint = 0, name = "log2OR",
                           limits = c(min(d$log2OR), max(d$log2OR)),
                           oob = squish) +
      scale_y_discrete(limits = rev(levels(d$cls)), labels = CLASS_LABELS) +
      labs(x = NULL, y = NULL, title = assay) +
      theme_pub() + theme(plot.margin = margin(4, 6, 12, 4),
                          axis.text.x = element_text(angle = 55, hjust = 1,
                                                     size = 4.8))
  }
  p11 <- tf_heat(me, "Motif") | tf_heat(be, "Binding")

  (p6 | p7) / (p8 | p10) / p11 +
    plot_annotation(tag_levels = "a", title = "C2  TE, conservation and TF") + theme(plot.margin = margin(4, 6, 12, 4))
}

# ---- C3: activity, chromatin and constraint --------------------------------
panel_c3 <- function() {
  ac <- read_res("cage_activity_by_class.tsv") %>%
    filter(cls %in% c("PCG", LNC5)) %>% mutate(cls = fct_cls(cls))
  p14 <- ggplot(ac, aes(cls, active_frac, fill = cls)) +
    geom_col(width = 0.72, colour = "black", linewidth = 0.15) +
    scale_fill_manual(values = CLASS_COLORS, guide = "none") +
    scale_x_discrete(labels = CLASS_LABELS) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = "CAGE-active", title = "Promoter activity") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))

  ch <- read_res("chromhmm_by_class.tsv") %>% filter(cls %in% c("PCG", LNC5))
  chm <- ch %>% select(cls, matches("_(prom_like|enh_like)$")) %>%
    pivot_longer(-cls, names_to = "column", values_to = "fraction") %>%
    mutate(cls = fct_cls(cls),
           state = if_else(str_detect(column, "prom_like"),
                           "Promoter-like", "Enhancer-like")) %>%
    group_by(cls, state) %>% summarise(fraction = mean(fraction), .groups = "drop")
  p15 <- ggplot(chm, aes(cls, fraction, fill = state)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7,
             colour = "black", linewidth = 0.15) +
    scale_fill_manual(values = c(`Promoter-like` = "#0072B2",
                                 `Enhancer-like` = "#E69F00"), name = NULL) +
    scale_x_discrete(labels = CLASS_LABELS) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Fraction of TSS", title = "ChromHMM state") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2),
                        legend.position = "top")

  cp <- read_res("constraint_structure", "A_conservation_tests.tsv") %>%
    filter(feature %in% c("prom_phylop", "body_phylop", "exon_phylop"),
           cls %in% c("PCG", LNC5)) %>%
    mutate(cls = fct_cls(cls),
           feature = factor(feature, levels = c("prom_phylop", "body_phylop",
                                                "exon_phylop"),
                            labels = c("Promoter", "Body", "Exon")),
           delta = median_class - median_PCG)
  p16 <- ggplot(cp, aes(cls, delta, fill = cls)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    facet_wrap(~feature, nrow = 1) +
    scale_fill_manual(values = CLASS_COLORS, guide = "none") +
    scale_x_discrete(labels = CLASS_LABELS) +
    labs(x = NULL, y = "phyloP delta", title = "Promoter-centred constraint") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))

  pt <- read_res("constraint_structure", "B_patch_tests.tsv") %>%
    filter(feature == "prom_patch_frac_phylop", cls %in% c("PCG", LNC5)) %>%
    mutate(cls = fct_cls(cls))
  p17 <- ggplot(pt, aes(cls, median_class, fill = cls)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.15) +
    geom_hline(yintercept = pt$median_PCG[1], linetype = 2, linewidth = 0.35) +
    scale_fill_manual(values = CLASS_COLORS, guide = "none") +
    scale_x_discrete(labels = CLASS_LABELS) +
    labs(x = NULL, y = "Patch fraction", title = "Conserved patches") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))

  td <- read_res("constraint_structure", "C_te_dual_tests.tsv") %>%
    filter(feature %in% c("te_cov", "te_exon_frac"), cls %in% c("PCG", LNC5)) %>%
    mutate(cls = fct_cls(cls),
           feature = factor(feature, levels = c("te_cov", "te_exon_frac"),
                            labels = c("Promoter TE", "TE exon")),
           delta = median_class - median_PCG)
  p18 <- ggplot(td, aes(cls, delta, fill = cls)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    facet_wrap(~feature, nrow = 1) +
    scale_fill_manual(values = CLASS_COLORS, guide = "none") +
    scale_x_discrete(labels = CLASS_LABELS) +
    labs(x = NULL, y = "Delta", title = "TE dual role") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4),
                        axis.text.x = element_text(angle = 90, hjust = 1,
                                                   vjust = 0.5, size = 4.0))

  tt <- read_res("ttseq", "tests.tsv") %>%
    filter(feature %in% c("log_syn", "log_decay", "log_ratio"),
           cls %in% c("PCG", LNC5)) %>%
    mutate(cls = fct_cls(cls),
           feature = factor(feature, levels = c("log_syn", "log_decay", "log_ratio"),
                            labels = c("Synthesis", "Decay", "Syn-decay")),
           delta = median_class - median_PCG)
  p19 <- ggplot(tt, aes(cls, delta, fill = cls)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    facet_wrap(~feature, nrow = 1) +
    scale_fill_manual(values = CLASS_COLORS, guide = "none") +
    scale_x_discrete(labels = CLASS_LABELS) +
    labs(x = NULL, y = "Delta", title = "TT-seq") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4),
                        axis.text.x = element_text(angle = 90, hjust = 1,
                                                   vjust = 0.5, size = 3.6))

  (p14 | p15) / (p16) / (p17 | p18 | p19) +
    plot_annotation(tag_levels = "a", title = "C3  Activity, chromatin and constraint") + theme(plot.margin = margin(4, 6, 12, 4))
}

# ---- C4: tissues, tumours and validation -----------------------------------
panel_c4 <- function() {
  gt <- read_res("tissue_disease", "tables", "gtex_class_by_tissue.tsv")
  # Only the five birth-class shares partition the expressed set; share_LTR is
  # an overlapping annotation and would double count inside the stack.
  share_cols <- paste0("share_", LNC5)
  gtl <- gt %>% select(tissue, all_of(share_cols), n_lnc_expressed) %>%
    pivot_longer(all_of(share_cols), names_to = "cls", values_to = "share") %>%
    mutate(cls = factor(str_remove(cls, "^share_"), levels = LNC5),
           tissue = fct_reorder(tissue, n_lnc_expressed))
  # 49 GTEx tissues cannot be labelled legibly in a composite-width panel, so
  # show the 20 with the most expressed lncRNAs and say so in the title.
  keep <- gtl %>% distinct(tissue, n_lnc_expressed) %>%
    slice_max(n_lnc_expressed, n = 20) %>% pull(tissue)
  gtl <- gtl %>% filter(tissue %in% keep) %>% mutate(tissue = fct_drop(tissue))
  p20 <- ggplot(gtl, aes(tissue, share, fill = cls)) +
    geom_col(width = 0.8, colour = "black", linewidth = 0.08) +
    coord_flip() +
    scale_fill_manual(values = CLASS_COLORS[LNC5],
                      labels = CLASS_LABELS[LNC5], name = NULL) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Share", title = "GTEx tissues (top 20)") +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4),
                        axis.text.y = element_text(size = 5),
                        legend.position = "bottom",
                        legend.text = element_text(size = 5.6),
                        legend.key.size = unit(2.5, "mm"))

  de <- read_res("tissue_disease", "tables", "tcga_class_de_per_cancer.tsv") %>%
    filter(cls %in% LNC5, !is.na(frac_up)) %>%
    mutate(cls = fct_lnc5(cls)) %>%
    group_by(cancer) %>% mutate(ord = mean(frac_up, na.rm = TRUE)) %>%
    ungroup() %>% mutate(cancer = fct_reorder(cancer, ord))
  p21 <- ggplot(de, aes(cancer, cls, fill = frac_up)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradient(low = "white", high = "#08306B", name = "Up") +
    scale_y_discrete(labels = CLASS_LABELS) +
    labs(x = NULL, y = NULL, title = "TCGA class DE") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 4.8))

  cg <- read_res("ccle_cgga", "tables", "cgga_class_by_group.tsv") %>%
    filter(group == "Grade", cls %in% c("LTR_derived_lnc", LNC5)) %>%
    mutate(cls = if_else(cls == "LTR_derived_lnc", "LTR_lnc", cls),
           level = factor(level, levels = c("WHO II", "WHO III", "WHO IV")))
  p22 <- ggplot(cg, aes(level, frac_expr, colour = cls, group = cls)) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.5) +
    scale_colour_manual(values = c(CLASS_COLORS[LNC5], LTR_lnc = "#000000"),
                        labels = c(CLASS_LABELS[LNC5], LTR_lnc = "LTR-lnc"),
                        name = NULL) +
    labs(x = NULL, y = "Detected", title = "CGGA grade") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(legend.position = "bottom",
                        legend.text = element_text(size = 4.8))

  cl <- read_res("ccle_cgga", "tables", "ccle_class_by_lineage.tsv") %>%
    filter(cls %in% LNC5) %>% mutate(cls = fct_lnc5(cls))
  top <- cl %>% distinct(lineage, n_lines) %>% slice_max(n_lines, n = 12) %>% pull(lineage)
  cl <- cl %>% filter(lineage %in% top)
  p23 <- ggplot(cl, aes(lineage, cls, fill = frac_expr)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradient(low = "white", high = "#08306B", name = "Frac") +
    scale_y_discrete(labels = CLASS_LABELS) +
    labs(x = NULL, y = NULL, title = "CCLE lineages") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 4.8))

  sc <- read_res("scrna", "class_by_celltype.tsv") %>%
    filter(cls %in% c("PCG", LNC5)) %>% mutate(cls = fct_cls(cls))
  p24 <- ggplot(sc, aes(cls, cell_type, fill = det_rate)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradient(low = "white", high = "#08306B", name = "Detect") +
    scale_x_discrete(labels = CLASS_LABELS) +
    labs(x = NULL, y = NULL, title = "scRNA cell types") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))

  mb <- read_res("cns_cohorts", "mb_per_sample.tsv") %>%
    select(sample, matches("_frac$")) %>%
    pivot_longer(-sample, names_to = "cls", values_to = "fraction") %>%
    mutate(cls = str_remove(cls, "_frac$")) %>% filter(cls %in% LNC5) %>%
    mutate(cls = fct_lnc5(cls))
  p25 <- ggplot(mb, aes(cls, fraction, fill = cls)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.25) +
    scale_fill_manual(values = CLASS_COLORS[LNC5], guide = "none") +
    scale_x_discrete(labels = CLASS_LABELS[LNC5]) +
    labs(x = NULL, y = "Fraction", title = "Medulloblastoma") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))

  sv <- read_res("evo_contexts", "A_survival_meta.tsv") %>% mutate(feature = fct_cls(feature))
  p26 <- ggplot(sv, aes(feature, z_meta, fill = z_meta > 0)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"),
                      guide = "none") +
    scale_x_discrete(labels = CLASS_LABELS) +
    labs(x = NULL, y = "Stouffer z", title = "Survival meta") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4),
                        axis.text.x = element_text(angle = 90, hjust = 1,
                                                   vjust = 0.5, size = 4.2))

  (p20 | p21) / (p22 | p23) / (p24 | p25 | p26) +
    plot_annotation(tag_levels = "a", title = "C4  Tissues, tumours and cohorts") + theme(plot.margin = margin(4, 6, 12, 4))
}

# ---- C5: evolution and strand parity ---------------------------------------
panel_c5 <- function() {
  age <- read_res("comparative", "E_class_x_age.tsv") %>%
    pivot_longer(-cls, names_to = "age", values_to = "n") %>%
    mutate(cls = fct_lnc5(cls),
           age = factor(age, levels = rev(AGE_LEVELS), labels = rev(AGE_LABELS)))
  p27 <- ggplot(age, aes(cls, n, fill = age)) +
    geom_col(width = 0.72, colour = "black", linewidth = 0.12) +
    scale_fill_manual(values = age_pal(), name = NULL) +
    scale_x_discrete(labels = CLASS_LABELS[LNC5]) +
    labs(x = NULL, y = "Loci", title = "Transcript age") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2),
                        legend.position = "top", legend.text = element_text(size = 4.8))

  tj <- read_res("evo_contexts", "B_tajima_tests.tsv") %>%
    filter(feature == "tajD_mean", cls %in% c("PCG", LNC5)) %>%
    mutate(cls = fct_cls(cls), delta = median_class - median_PCG)
  p28 <- ggplot(tj, aes(cls, delta, fill = cls)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    scale_fill_manual(values = CLASS_COLORS, guide = "none") +
    scale_x_discrete(labels = CLASS_LABELS) +
    labs(x = NULL, y = "Tajima's D delta", title = "Recent selection") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))

  lf <- read_res("ltr_family", "family_by_class.tsv") %>%
    filter(cls %in% LNC5) %>% mutate(cls = fct_lnc5(cls))
  p29 <- ggplot(lf, aes(family, cls, fill = log2OR)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                         midpoint = 0, name = "log2OR",
                         na.value = "grey90") +
    scale_y_discrete(labels = CLASS_LABELS[LNC5]) +
    labs(x = NULL, y = NULL, title = "LTR families (grey, not tested)") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4.8))

  mc <- read_res("mouse_taxonomy", "cross_species_class_comparison.tsv") %>%
    pivot_longer(c(human_fraction, mouse_fraction), names_to = "species",
                 values_to = "fraction") %>%
    mutate(cls = fct_cls(cls),
           species = factor(species, levels = c("human_fraction", "mouse_fraction"),
                            labels = c("Human", "Mouse")))
  p30 <- ggplot(mc, aes(cls, fraction, fill = species)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7,
             colour = "black", linewidth = 0.15) +
    scale_fill_manual(values = c(Human = "#0072B2", Mouse = "#E69F00"), name = NULL) +
    scale_x_discrete(labels = CLASS_LABELS) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Fraction", title = "Mouse vs human") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2),
                        legend.position = "top")

  cs <- read_res("developmental", "cross_species.tsv") %>%
    filter(cls %in% c("PCG", LNC5)) %>%
    mutate(cls = fct_cls(cls), species = factor(species, levels = SPECIES_ORDER))
  # Only intergenic lncRNAs reach the syntenic breadth test outside mouse, so
  # n is printed only for that row; PCG cells are too narrow for a long label.
  # rho spans 0.73-0.86, so the fill is sequential, not diverging at zero.
  p31 <- ggplot(cs, aes(species, cls, fill = rho_breadth)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = if_else(cls == "intergenic",
                                  sprintf("%.2f (n=%d)", rho_breadth, n),
                                  sprintf("%.2f", rho_breadth))),
              size = 1.6) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08306B",
                        limits = c(0.70, 0.87), name = "rho") +
    scale_y_discrete(labels = c(PCG = "PCG (reference)",
                                intergenic = "Intergenic")) +
    labs(x = NULL, y = NULL, title = "Breadth conservation") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4),
                        axis.text.x = element_text(angle = 45, hjust = 1,
                                                   size = 4.8))

  ss <- read_res("strand_controls", "S_sources.tsv") %>%
    mutate(source = fct_reorder(source, plus_frac))
  p32 <- ggplot(ss, aes(source, plus_frac)) +
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, linewidth = 0.35) +
    geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.35) +
    coord_flip() + scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Plus fraction", title = "Strand balance") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.y = element_text(size = 4.8))

  lp <- read_res("strand_parity", "C_length_parity.tsv") %>%
    filter(feature == "gene", cls %in% c("PCG", LNC5)) %>% mutate(cls = fct_cls(cls))
  p33 <- ggplot(lp, aes(cls, ks_stat, fill = cls)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.15) +
    scale_fill_manual(values = CLASS_COLORS, guide = "none") +
    scale_x_discrete(labels = CLASS_LABELS) +
    labs(x = NULL, y = "KS statistic", title = "Length parity") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))

  (p27 | p28) / (p29 | p30) / (p31 | p32 | p33) +
    plot_annotation(tag_levels = "a", title = "C5  Evolution and strand parity") + theme(plot.margin = margin(4, 6, 12, 4))
}

# ---- C6: ML and perturbation -----------------------------------------------
panel_c6 <- function() {
  cm <- suppressMessages(read_res("ml", "A_confusion.tsv"))
  names(cm)[1] <- "true_class"   # first column is unnamed in the TSV
  cm <- cm %>%
    pivot_longer(-true_class, names_to = "predicted", values_to = "n") %>%
    group_by(true_class) %>% mutate(fraction = n / sum(n)) %>% ungroup() %>%
    mutate(true_class = fct_lnc5(true_class), predicted = fct_lnc5(predicted))
  p34 <- ggplot(cm, aes(predicted, true_class, fill = fraction)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.2f", fraction)), size = 2.8,
              colour = if_else(cm$fraction > 0.5, "white", "black")) +
    scale_fill_gradient(low = "white", high = "#08306B", name = "Frac") +
    scale_x_discrete(labels = CLASS_LABELS[LNC5]) +
    scale_y_discrete(limits = rev, labels = CLASS_LABELS[LNC5]) +
    labs(x = "Predicted", y = "True", title = "Class prediction") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))

  ab <- read_res("evo_ext", "A_ablation.tsv")
  p35 <- ggplot(ab, aes(feature_set, macro_f1, fill = delta_vs_F0 > 0)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.15) +
    geom_text(aes(label = sprintf("%.3f", macro_f1)), vjust = -0.4, size = 2.2) +
    scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#BBBBBB"),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Macro-F1", title = "phyloP ablation") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2))

  im <- read_res("evo_ext", "A_feature_importance_F2.tsv") %>%
    slice_max(abs(delta_f1), n = 12) %>% mutate(feature = fct_reorder(feature, delta_f1))
  p36 <- ggplot(im, aes(feature, delta_f1, fill = delta_f1 > 0)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.15) +
    coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#BBBBBB"),
                      guide = "none") +
    labs(x = NULL, y = "Delta F1", title = "Feature importance") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4))

  pt <- read_res("perturb", "perturb_tests.tsv") %>%
    filter(contrast == "LTR_lnc_vs_other") %>%
    mutate(dataset_orig = dataset,
           dataset_label = factor(dataset,
             levels = c("GSE85461_iDUX4_fastq", "GSE85461_lentiDUX4_14h",
                        "GSE56787_FSHD_vs_ctrl", "GSE337536_decitabine",
                        "GSE305426_DNMTi_vs_DMSO"),
             labels = c("iDUX4", "lenti-DUX4", "FSHD", "Decitabine", "DAC/AZA")),
           direction = case_when(str_detect(dataset_orig, "GSE85461") ~ "DUX4",
                                 str_detect(dataset_orig, "FSHD") ~ "Endogenous",
                                 TRUE ~ "DNMTi"))
  p37 <- ggplot(pt, aes(dataset_label, cliffs, fill = direction)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    geom_text(aes(label = star_lab(pval),
                  vjust = if_else(cliffs > 0, -0.4, 1.1)), size = 2.4) +
    scale_fill_manual(values = c(DUX4 = "#D55E00", Endogenous = "#999999",
                                 DNMTi = "#0072B2"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) +
    labs(x = NULL, y = "Cliff's delta", title = "Perturbation validation") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4)) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.2),
                        legend.position = "top")

  ps <- read_res("perturbseq", "class_tests_within_target.tsv") %>%
    filter(contrast == "divergent_vs_other_lnc") %>%
    mutate(target = fct_reorder(target, cliffs))
  p38 <- ggplot(ps, aes(target, cliffs, fill = qval < 0.05)) +
    geom_col(width = 0.6, colour = "black", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#BBBBBB"),
                      guide = "none") +
    labs(x = NULL, y = "Cliff's delta", title = "Perturb-seq") +
    theme_pub() + theme(plot.margin = margin(4, 6, 12, 4))

  (p34 | p35) / (p36 | p37) / p38 +
    plot_annotation(tag_levels = "a", title = "C6  Machine learning and perturbation") + theme(plot.margin = margin(4, 6, 12, 4))
}

# ---- build and save ---------------------------------------------------------
c1 <- panel_c1(); save_composite(c1, "C1_classes_promoters", 183, 230)
c2 <- panel_c2(); save_composite(c2, "C2_te_conservation_tf", 183, 250)
c3 <- panel_c3(); save_composite(c3, "C3_activity_constraint", 183, 240)
c4 <- panel_c4(); save_composite(c4, "C4_tissues_tumours", 183, 250)
c5 <- panel_c5(); save_composite(c5, "C5_evolution_parity", 183, 240)
c6 <- panel_c6(); save_composite(c6, "C6_ml_perturbation", 183, 220)

message("Composite figures written to ", FIGDIR)
