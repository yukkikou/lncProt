# 15 final gap-filling statistics

_Generated 2026-10-02 09:45:02 by 15_final_gaps_

- FANTOM~GTEx tau: rho=nan, p=2.9e-137, n=23746
- Testis LTR enrichment: OR=0.27, p=2.9e-137 (per-tissue aggregate)
- Pan-cancer LTR derepression: 0.005 [-0.025,0.036], P=0.383

## TF binding vs LTR-derived promoters

| tf | bound_lnc | bound_LTR_derived | frac_LTR_among_bound | background_LTR_frac | fisher_OR | p |
|---|---|---|---|---|---|---|
| DUX4 | 278 | 103 | 0.371 | 0.218 | 2.149 | 0.000 |
| PHOX2B | 522 | 55 | 0.105 | 0.218 | 0.411 | 0.000 |
| BATF3 | 159 | 24 | 0.151 | 0.218 | 0.634 | 0.042 |
| CTCF | 7036 | 993 | 0.141 | 0.218 | 0.389 | 0.000 |
| MAZ | 4031 | 362 | 0.090 | 0.218 | 0.266 | 0.000 |

## Class activity vs dysregulation

| cls | gtex_median_frac_expr | tcga_median_frac_up |
|---|---|---|
| divergent | 0.415 | 0.121 |
| antisense | 0.354 | 0.135 |
| intronic | 0.193 | 0.098 |
| sense_overlapping | 0.371 | 0.117 |
| intergenic | 0.146 | 0.095 |
rho=0.70 (p=0.19)
