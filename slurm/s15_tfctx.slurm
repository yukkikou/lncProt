#!/bin/bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH -J lncP-tfctx
#SBATCH -p cpu
#SBATCH -D /media/share/workdir/1610305236/lncProt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH -o /media/share/home/1610305236/log/%x-%j.out
#SBATCH -e /media/share/home/1610305236/log/%x-%j.err
source /media/share/home/1610305236/.local/share/mamba/etc/profile.d/mamba.sh
micromamba activate Jupyter
run() { local m=$1; shift
    if [ -f "$m" ]; then echo ">>> skip: $m"; return 0; fi
    "$@" && touch "$m" || { echo "FAILED: $*" >&2; exit 1; } }
run data/.done_31t python -W ignore scripts/31_tf_contexts.py T
echo "s15 done: $(date)"
