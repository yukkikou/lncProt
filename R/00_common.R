# =============================================================================
# 00_common.R  --  lncProt unified plotting configuration
# -----------------------------------------------------------------------------
# Run from the project root (or set LNCPROT_DIR):
#   Sys.setenv(LNCPROT_DIR = "/path/to/lncProt")
#   source("R/00_common.R")
#
# Everything downstream uses these shared definitions:
#   * CLASS_ORDER / CLASS_COLORS / CLASS_LABELS  -- all birth classes + PCG
#   * LNC5                                        -- the five lncRNA classes
#   * AGE_LEVELS / AGE_COLORS                     -- transcript-age strata
#   * TE_COLORS / SPECIES_COLORS                  -- auxiliary palettes
#   * theme_pub()                                 -- publication theme
#   * save_panel()                                -- PDF + 600 dpi PNG
#   * read_res()                                  -- read results/<...>.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

PROJ <- Sys.getenv("LNCPROT_DIR", unset = normalizePath("."))
if (!dir.exists(file.path(PROJ, "results"))) {
  stop("Cannot find results/ under PROJ = ", PROJ,
       ". Run from the lncProt project root or set LNCPROT_DIR.")
}
RESDIR <- file.path(PROJ, "results")
FIGDIR <- file.path(PROJ, "figures", "R")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

read_res <- function(...) {
  readr::read_tsv(file.path(RESDIR, ...), show_col_types = FALSE)
}

# ---- unified group definitions ---------------------------------------------
CLASS_ORDER <- c("PCG", "divergent", "antisense", "intronic",
                 "sense_overlapping", "intergenic", "pseudogene",
                 "transcribed_pseudogene")
LNC5 <- c("divergent", "antisense", "intronic", "sense_overlapping",
          "intergenic")

# Okabe-Ito, colour-blind safe, with two substitutions checked under
# deuteranopia. Orange and vermillion (original divergent and intronic) sit at
# dE 12.8 and were too close, so intronic is violet and pseudogene a darker grey.
# Keep identical across every script.
CLASS_COLORS <- c(
  PCG                 = "#0072B2",
  divergent           = "#E69F00",
  antisense           = "#009E73",
  intronic            = "#5D3A9B",
  sense_overlapping   = "#CC79A7",
  intergenic          = "#56B4E9",
  pseudogene          = "#555555",
  transcribed_pseudogene = "#333333"
)

CLASS_LABELS <- c(
  PCG = "PCG", divergent = "Divergent", antisense = "Antisense",
  intronic = "Intronic", sense_overlapping = "Sense-overlap",
  intergenic = "Intergenic", pseudogene = "Pseudogene",
  transcribed_pseudogene = "Transcribed pseudo",
  LTR_lnc = "LTR-lnc", other_lnc = "Other lnc"
)

AGE_LEVELS <- c("human_specific", "primate", "mammalian", "eutherian")
AGE_LABELS <- c(human_specific = "Human-specific", primate = "Primate",
                mammalian = "Mammalian", eutherian = "Eutherian")
AGE_COLORS <- c(human_specific = "#D9D9D9", primate = "#56B4E9",
                mammalian = "#0072B2", eutherian = "#003B5C")
# Display labels are what the factor levels carry in plots, so the palette must
# be re-keyed by label before it is handed to scale_fill_manual().
age_pal <- function() setNames(unname(AGE_COLORS[AGE_LEVELS]),
                               unname(AGE_LABELS[AGE_LEVELS]))

TE_LEVELS <- c("LTR", "LINE", "SINE", "DNA")
TE_COLORS <- c(LTR = "#D55E00", LINE = "#0072B2",
               SINE = "#009E73", DNA = "#CC79A7")

SPECIES_ORDER <- c("human", "macaque", "mouse", "rat", "rabbit",
                   "chicken", "opossum")
SPECIES_COLORS <- c(human = "#000000", macaque = "#E69F00",
                    mouse = "#0072B2", rat = "#56B4E9",
                    rabbit = "#009E73", chicken = "#CC79A7",
                    opossum = "#999999")

# ordered class factor helper; unknown values are kept but sorted last
fct_cls <- function(x) {
  factor(x, levels = c(CLASS_ORDER, setdiff(unique(x), CLASS_ORDER)))
}
fct_lnc5 <- function(x) factor(x, levels = LNC5)

# ---- publication theme ------------------------------------------------------
theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(linewidth = 0.35),
      axis.line = element_line(linewidth = 0.4),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size - 0.5),
      plot.tag = element_text(face = "bold", size = base_size + 1),
      plot.title = element_text(size = base_size + 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(4, 6, 4, 4)
    )
}

# Save both vector PDF and 600-dpi PNG. width/height in mm.
save_panel <- function(plot, name, width = 180, height = 120, dpi = 600) {
  ggsave(file.path(FIGDIR, paste0(name, ".pdf")), plot,
         width = width / 25.4, height = height / 25.4, units = "in")
  ggsave(file.path(FIGDIR, paste0(name, ".png")), plot,
         width = width / 25.4, height = height / 25.4, units = "in",
         dpi = dpi)
  invisible(plot)
}

# common significance labels used in several panels
star_lab <- function(fdr) {
  case_when(
    is.na(fdr) ~ "",
    fdr < 1e-10 ~ "***",
    fdr < 1e-3  ~ "**",
    fdr < 0.05  ~ "*",
    TRUE        ~ "ns"
  )
}
