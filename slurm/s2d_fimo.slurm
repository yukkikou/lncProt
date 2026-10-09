#!/bin/bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH -J lncP-fimo
#SBATCH -p cpu
#SBATCH -D /media/share/workdir/1610305236/lncProt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --array=0-9%5
#SBATCH -o /media/share/home/1610305236/log/%x-%A-%a.out
#SBATCH -e /media/share/home/1610305236/log/%x-%A-%a.err

source /media/share/home/1610305236/.local/share/mamba/etc/profile.d/mamba.sh
micromamba activate Jupyter

CHUNK=data/motifs/chunks/chunk_${SLURM_ARRAY_TASK_ID}.meme
OK=data/.done_fimo_${SLURM_ARRAY_TASK_ID}
if [ -f "$OK" ]; then echo ">>> skip fimo chunk $SLURM_ARRAY_TASK_ID"; exit 0; fi
if [ "$SLURM_ARRAY_TASK_ID" = "0" ] && [ ! -f data/motifs/chunks/chunk_9.meme ]; then
    bash scripts/07a_fimo_prepare.sh prepare || { echo "FIMO prepare failed" >&2; exit 1; }
fi
sleep 10
[ -f "$CHUNK" ] || { echo "missing $CHUNK" >&2; exit 1; }
bash scripts/07a_fimo_prepare.sh scan && touch "$OK" || { echo "FIMO scan chunk failed" >&2; exit 1; }
echo "s2d chunk $SLURM_ARRAY_TASK_ID complete: $(date)"
