#!/bin/bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH -J lncP-td
#SBATCH -p cpu
#SBATCH -D /media/share/workdir/1610305236/lncProt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH -o /media/share/home/1610305236/log/%x-%j.out
#SBATCH -e /media/share/home/1610305236/log/%x-%j.err

# s5 driver: parse Toil -> GTEx tissues -> (array job handles per-cancer TCGA)
# -> merge -> stage -> integrate
source /media/share/home/1610305236/.local/share/mamba/etc/profile.d/mamba.sh
micromamba activate Jupyter

run() { local marker=$1; shift; if [ -f "$marker" ]; then echo ">>> skip $marker"; return 0; fi
    "$@" && touch "$marker" || { echo "FAILED: $*" >&2; exit 1; }; }

# 13b: parse + accounting (matrix already downloaded by 13a)
run data/tissue_disease/.done_13b python -W ignore scripts/13b_parse_toil.py

# 13c: GTEx tissue dimension
run data/tissue_disease/.done_13c python -W ignore scripts/13c_gtex_tissue.py

# 13d per-cancer: submit the ARRAY from here and wait, so the whole chain is
# one SLURM dependency graph
CANCERS=$(python - <<'PYEOF'
import pandas as pd
s = pd.read_csv("data/tissue_disease/toil_samples.tsv", sep="\t")
t = s[(s.cohort=="TCGA") & (s.condition=="tumor")].groupby("cancer").size()
print(" ".join(t[t>=20].index))
PYEOF
)
echo "cancers with n>=20 tumors: $CANCERS"
ARR=$(sbatch --parsable --array=0-$(($(wc -w <<< "$CANCERS")-1))%8 \
      -J lncP-tdc --mem=24G --time=04:00:00 \
      -o /media/share/home/1610305236/log/lncP-tdc-%A-%a.out \
      -e /media/share/home/1610305236/log/lncP-tdc-%A-%a.err \
      --wrap "source /media/share/home/1610305236/.local/share/mamba/etc/profile.d/mamba.sh; micromamba activate Jupyter; cd $PWD; C=($CANCERS); python -W ignore scripts/13d_tcga_disease.py one \${C[\$SLURM_ARRAY_TASK_ID]}")
echo "array submitted: $ARR"
sbatch --wait $ARR || true

# merge + stage + integrate
run data/tissue_disease/.done_13dm python -W ignore scripts/13d_tcga_disease.py merge
run data/tissue_disease/.done_13ds python -W ignore scripts/13d_tcga_disease.py stage
run data/tissue_disease/.done_13e python -W ignore scripts/13e_integrate_td.py
run data/tissue_disease/.done_14  python -W ignore scripts/14_td_figures.py

echo "s5 complete: $(date)"
