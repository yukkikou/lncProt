# 12d robustness controls

_Generated 2026-10-01 10:08:41 by 12d_robustness_

## Divergent window sensitivity

| 500bp | 1000bp (primary) | 2000bp |
|---|---|---|
| 2118 | 2692 | 3201 |

## Active-promoter-only comparison

| feature | cls | n_active_class | mean_active_class | mean_active_PCG | diff |
|---|---|---|---|---|---|
| phastCons | divergent | 1668 | 0.2844 | 0.2550 | 0.0294 |
| phastCons | antisense | 783 | 0.2066 | 0.2550 | -0.0484 |
| phastCons | intronic | 695 | 0.1657 | 0.2550 | -0.0893 |
| phastCons | intergenic | 3270 | 0.1283 | 0.2550 | -0.1267 |
| prom-like chromatin frac | divergent | 1668 | 0.7857 | 0.5924 | 0.1933 |
| prom-like chromatin frac | antisense | 783 | 0.2450 | 0.5924 | -0.3475 |
| prom-like chromatin frac | intronic | 695 | 0.1301 | 0.5924 | -0.4624 |
| prom-like chromatin frac | intergenic | 3271 | 0.1339 | 0.5924 | -0.4585 |
| GC content | divergent | 1668 | 0.6470 | 0.6108 | 0.0362 |
| GC content | antisense | 783 | 0.5412 | 0.6108 | -0.0695 |
| GC content | intronic | 695 | 0.4836 | 0.6108 | -0.1272 |
| GC content | intergenic | 3271 | 0.5028 | 0.6108 | -0.1080 |

## Conservation vs class, GC-controlled (partial Spearman)

| cls | n | rho_raw | p_raw | rho_GC_controlled | p_controlled |
|---|---|---|---|---|---|
| divergent | 22757 | 0.046 | 0.000 | -0.003 | 0.693 |
| antisense | 22841 | -0.140 | 0.000 | -0.033 | 0.000 |
| intronic | 23464 | -0.283 | 0.000 | -0.124 | 0.000 |
| sense_overlapping | 20415 | -0.056 | 0.000 | 0.011 | 0.116 |
| intergenic | 30089 | -0.424 | 0.000 | -0.283 | 0.000 |
