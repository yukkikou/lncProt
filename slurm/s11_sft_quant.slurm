#!/bin/bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH -J lncP-sftq
#SBATCH -p cpu
#SBATCH -D /media/share/workdir/1610305236/lncProt
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=0-106%16
#SBATCH -o /media/share/home/1610305236/log/lncP-sftq-%A-%a.out
#SBATCH -e /media/share/home/1610305236/log/lncP-sftq-%A-%a.err
W=/media/share/workdir/1610305236/lncProt
SRC=/media/share/data8/volume1/1610305236_CJ_SFT/RNA-seq
S=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $W/data/sft_own/sample_ids.txt)
OUT=$W/data/sft_own/quants/$S
[ -s $OUT/quant.sf ] && { echo "skip $S"; exit 0; }
mkdir -p $OUT
singularity exec --no-mount bind-paths --bind /media/share/data8/volume1:/media/share/data8/volume1:ro --bind $W:$W \
  /media/share/resource/container/salmon/salmon-1.10.2.sif salmon quant \
  -i $W/data/gbm_own/ref/salmon_idx -l A -p 8 \
  -1 $SRC/${S}_1.fq.gz -2 $SRC/${S}_2.fq.gz \
  --validateMappings -o $OUT
echo "$S done: $(date)"
