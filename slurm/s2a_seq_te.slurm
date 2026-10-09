#!/bin/bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH -J lncP-seqte
#SBATCH -p cpu
#SBATCH -D /media/share/workdir/1610305236/lncProt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH -o /media/share/home/1610305236/log/%x-%j.out
#SBATCH -e /media/share/home/1610305236/log/%x-%j.err

source /media/share/home/1610305236/.local/share/mamba/etc/profile.d/mamba.sh
micromamba activate Jupyter

run() { # run <marker> <cmd...> — restartable + fail-fast
    local marker=$1; shift
    if [ -f "$marker" ]; then echo ">>> skip (done): $marker"; return 0; fi
    "$@" && touch "$marker" || { echo "FAILED: $*" >&2; exit 1; }
}

run data/.done_05 python -W ignore scripts/05_seq_features.py
run data/.done_06 python -W ignore scripts/06_te_promoters.py
echo "s2a complete: $(date)"
