#!/bin/bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH -J lncP-phylop
#SBATCH -p cpu
#SBATCH -D /media/share/workdir/1610305236/lncProt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH -o /media/share/home/1610305236/log/%x-%j.out
#SBATCH -e /media/share/home/1610305236/log/%x-%j.err

source /media/share/home/1610305236/.local/share/mamba/etc/profile.d/mamba.sh
micromamba activate Jupyter

run() { # run <marker> <cmd...> — restartable + fail-fast
    local marker=$1; shift
    if [ -f "$marker" ]; then echo ">>> skip (done): $marker"; return 0; fi
    "$@" && touch "$marker" || { echo "FAILED: $*" >&2; exit 1; }
}

run data/.done_27 python -W ignore scripts/27_phylop.py
# figure legs: uncomment after 16/11 fig6 code is edited for phyloP panels
run data/.done_16b python -W ignore scripts/16_submission_figures.py
run data/.done_11b python -W ignore scripts/11_report_figures.py
echo "s12 phylop complete: $(date)"
