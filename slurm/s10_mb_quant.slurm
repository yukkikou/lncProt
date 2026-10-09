#!/bin/bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH -J lncP-mbq
#SBATCH -p cpu
#SBATCH -D /media/share/workdir/1610305236/lncProt
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=0-205%16
#SBATCH -o /media/share/home/1610305236/log/lncP-mbq-%A-%a.out
#SBATCH -e /media/share/home/1610305236/log/lncP-mbq-%A-%a.err
W=/media/share/workdir/1610305236/lncProt
SRC_DIR=$(grep -l "" /dev/null; dirname "$(find /media/share/data8/volume1/1610305236_CJ_MB -name "$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $W/data/mb_own/sample_ids.txt)_R1.fastq.gz" | head -1)")
S=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $W/data/mb_own/sample_ids.txt)
OUT=$W/data/mb_own/quants/$S
[ -s $OUT/quant.sf ] && { echo "skip $S"; exit 0; }
mkdir -p $OUT
singularity exec --no-mount bind-paths --bind /media/share/data8/volume1:/media/share/data8/volume1:ro --bind $W:$W \
  /media/share/resource/container/salmon/salmon-1.10.2.sif salmon quant \
  -i $W/data/gbm_own/ref/salmon_idx -l A -p 8 \
  -1 $SRC_DIR/${S}_R1.fastq.gz -2 $SRC_DIR/${S}_R2.fastq.gz \
  --validateMappings -o $OUT
echo "$S done: $(date)"
