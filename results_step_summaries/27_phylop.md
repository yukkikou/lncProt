# 27 phyloP 241-mammal conservation + acceleration

_Generated 2026-10-08 15:20:28 by 27_phylop_

## A. Class vs PCG (phyloP metrics)

| feature | cls | n_class | mean_class | mean_PCG | cliffs_delta | pval | fdr |
|---|---|---|---|---|---|---|---|
| frac accelerated (phylop<=-2) | antisense | 2776 | 0.0483 | 0.0578 | -0.1862 | 0.0000 | 0.0000 |
| frac accelerated (phylop<=-2) | divergent | 2691 | 0.0608 | 0.0578 | 0.0400 | 0.0007 | 0.0007 |
| frac accelerated (phylop<=-2) | intergenic | 10006 | 0.0443 | 0.0578 | -0.2613 | 0.0000 | 0.0000 |
| frac accelerated (phylop<=-2) | intronic | 3401 | 0.0435 | 0.0578 | -0.2745 | 0.0000 | 0.0000 |
| frac accelerated (phylop<=-2) | pseudogene | 13549 | 0.0322 | 0.0578 | -0.4828 | 0.0000 | 0.0000 |
| frac accelerated (phylop<=-2) | sense_overlapping | 350 | 0.0500 | 0.0578 | -0.1406 | 0.0000 | 0.0000 |
| frac accelerated (phylop<=-2) | transcribed_pseudogene | 1630 | 0.0446 | 0.0578 | -0.2468 | 0.0000 | 0.0000 |
| frac conserved (phylop>=+2) | antisense | 2776 | 0.1362 | 0.1979 | -0.3265 | 0.0000 | 0.0000 |
| frac conserved (phylop>=+2) | divergent | 2691 | 0.2264 | 0.1979 | 0.0874 | 0.0000 | 0.0000 |
| frac conserved (phylop>=+2) | intergenic | 10006 | 0.0673 | 0.1979 | -0.6301 | 0.0000 | 0.0000 |
| frac conserved (phylop>=+2) | intronic | 3401 | 0.0749 | 0.1979 | -0.5896 | 0.0000 | 0.0000 |
| frac conserved (phylop>=+2) | pseudogene | 13549 | 0.0611 | 0.1979 | -0.6395 | 0.0000 | 0.0000 |
| frac conserved (phylop>=+2) | sense_overlapping | 350 | 0.1363 | 0.1979 | -0.3033 | 0.0000 | 0.0000 |
| frac conserved (phylop>=+2) | transcribed_pseudogene | 1630 | 0.1055 | 0.1979 | -0.4087 | 0.0000 | 0.0000 |
| phyloP (241w) mean | antisense | 2776 | 0.5578 | 0.7789 | -0.2065 | 0.0000 | 0.0000 |
| phyloP (241w) mean | divergent | 2691 | 0.9137 | 0.7789 | 0.0576 | 0.0000 | 0.0000 |
| phyloP (241w) mean | intergenic | 10006 | 0.1682 | 0.7789 | -0.4827 | 0.0000 | 0.0000 |
| phyloP (241w) mean | intronic | 3401 | 0.2325 | 0.7789 | -0.4224 | 0.0000 | 0.0000 |
| phyloP (241w) mean | pseudogene | 13549 | 0.2426 | 0.7789 | -0.3989 | 0.0000 | 0.0000 |
| phyloP (241w) mean | sense_overlapping | 350 | 0.5254 | 0.7789 | -0.2133 | 0.0000 | 0.0000 |
| phyloP (241w) mean | transcribed_pseudogene | 1630 | 0.3972 | 0.7789 | -0.2580 | 0.0000 | 0.0000 |

## B. GC-controlled partial Spearman (phylop_mean)

| cls | n | rho_raw | p_raw | rho_GC_controlled | p_controlled |
|---|---|---|---|---|---|
| divergent | 22715 | 0.032 | 0.000 | -0.001 | 0.893 |
| antisense | 22800 | -0.117 | 0.000 | -0.027 | 0.000 |
| intronic | 23425 | -0.258 | 0.000 | -0.121 | 0.000 |
| sense_overlapping | 20374 | -0.048 | 0.000 | 0.006 | 0.358 |
| intergenic | 30030 | -0.394 | 0.000 | -0.278 | 0.000 |

## C. LTR-origin acceleration contrasts

| contrast | feature | n_case | n_ctrl | median_case | median_ctrl | mean_case | mean_ctrl | cliffs_delta | pval | fdr |
|---|---|---|---|---|---|---|---|---|---|---|
| lnc LTR-TSS vs non-LTR | frac_accel_2 | 1584 | 17640 | 0.0243 | 0.0400 | 0.0319 | 0.0485 | -0.3184 | 0.0000 | 0.0000 |
| lnc LTR-TSS vs non-LTR | phylop_mean | 1584 | 17640 | -0.0212 | 0.0771 | -0.0103 | 0.3787 | -0.2687 | 0.0000 | 0.0000 |
| intergenic LTR-TSS vs non-LTR | frac_accel_2 | 1172 | 8834 | 0.0229 | 0.0371 | 0.0314 | 0.0461 | -0.2881 | 0.0000 | 0.0000 |
| intergenic LTR-TSS vs non-LTR | phylop_mean | 1172 | 8834 | -0.0293 | 0.0082 | -0.0314 | 0.1946 | -0.1496 | 0.0000 | 0.0000 |

## D. phastCons vs phyloP by class

| cls | n | phast_mean | phylop_mean | frac_accel_2 | frac_cons_2 | n_ltr | frac_accel2_ltr | spearman_phast_phylop | p_spearman |
|---|---|---|---|---|---|---|---|---|---|
| PCG | 20065 | 0.242 | 0.779 | 0.058 | 0.198 | 469 | 0.054 | 0.931 | 0.000 |
| divergent | 2692 | 0.274 | 0.914 | 0.061 | 0.226 | 10 |  | 0.960 | 0.000 |
| antisense | 2778 | 0.185 | 0.558 | 0.048 | 0.136 | 153 | 0.033 | 0.909 | 0.000 |
| intronic | 3403 | 0.120 | 0.232 | 0.044 | 0.075 | 223 | 0.033 | 0.839 | 0.000 |
| sense_overlapping | 350 | 0.178 | 0.525 | 0.050 | 0.136 | 28 | 0.030 | 0.888 | 0.000 |
| intergenic | 10035 | 0.108 | 0.168 | 0.044 | 0.067 | 1181 | 0.031 | 0.772 | 0.000 |
| pseudogene | 13577 | 0.174 | 0.243 | 0.032 | 0.061 | 1373 | 0.028 | 0.449 | 0.000 |
| transcribed_pseudogene | 1633 | 0.166 | 0.397 | 0.045 | 0.106 | 109 | 0.033 | 0.620 | 0.000 |
| all_promoters | 54533 | 0.189 | 0.481 | 0.047 | 0.127 | 3546 | 0.033 | 0.805 | 0.000 |
