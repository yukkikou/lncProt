#!/bin/bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH -J lncP-cc
#SBATCH -p cpu
#SBATCH -D /media/share/workdir/1610305236/lncProt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH -o /media/share/home/1610305236/log/%x-%j.out
#SBATCH -e /media/share/home/1610305236/log/%x-%j.err

# s6: CCLE/CGGA extension — runs when matrices land in data/ccle_cgga/
source /media/share/home/1610305236/.local/share/mamba/etc/profile.d/mamba.sh
micromamba activate Jupyter

if ls data/ccle_cgga/ccle/*.csv data/ccle_cgga/cgga/*.txt* >/dev/null 2>&1; then
    python -W ignore scripts/19_ccle_cgga.py
else
    echo "matrices not yet present — see data/ccle_cgga/README_PLACE_FILES_HERE.md"
    exit 1
fi
