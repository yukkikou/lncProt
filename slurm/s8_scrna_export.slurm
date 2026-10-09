#!/bin/bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH -J lncP-scrna
#SBATCH -p cpu
#SBATCH -D /media/share/workdir/1610305236/lncProt
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH -o /media/share/home/1610305236/log/%x-%j.out
#SBATCH -e /media/share/home/1610305236/log/%x-%j.err
REF=/media/share/resource/reference
W=/media/share/workdir/1610305236/lncProt
singularity exec --no-mount bind-paths --bind $REF:$REF:ro --bind $W:$W \
  /media/share/resource/container/R_env/seurat/seurat.sif Rscript $W/scripts/21_scrna_export.R
