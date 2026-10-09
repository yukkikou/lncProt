#!/bin/bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH -J lncP-prep
#SBATCH -p cpu
#SBATCH -D /media/share/workdir/1610305236/lncProt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH -o /media/share/home/1610305236/log/%x-%j.out
#SBATCH -e /media/share/home/1610305236/log/%x-%j.err

# s1: reference manifest + GENCODE parse + ChromHMM liftover + locus classification
source /media/share/home/1610305236/.local/share/mamba/etc/profile.d/mamba.sh
micromamba activate Jupyter

run() { # run <marker> <cmd...> — restartable + fail-fast
    local marker=$1; shift
    if [ -f "$marker" ]; then echo ">>> skip (done): $marker"; return 0; fi
    "$@" && touch "$marker" || { echo "FAILED: $*" >&2; exit 1; }
}

source scripts/common.sh
progress_init
ref_manifest build && echo "manifest built"

run data/.done_02 python -W ignore scripts/02_parse_gencode.py

echo ">>> Step 03b: ChromHMM liftover + QC"
[ -f data/chromhmm_hg38/.done ] || bash scripts/03b_chromhmm.sh || echo "WARN: 03b failed - chromatin module degraded"

run data/.done_04 python -W ignore scripts/04_classify_loci.py
run data/.done_04b python -W ignore scripts/04b_lnc_antisense_flag.py

echo "s1 complete: $(date)"
