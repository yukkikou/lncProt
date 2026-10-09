#!/bin/bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH -J lncP-idx
#SBATCH -p cpu
#SBATCH -D /media/share/workdir/1610305236/lncProt
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH -o /media/share/home/1610305236/log/%x-%j.out
#SBATCH -e /media/share/home/1610305236/log/%x-%j.err
set -e
REF=/media/share/resource/reference
GTF=$REF/homo_sapiens/GRCh38/annotation/gencode.v46.annotation.gtf
GEN=$REF/homo_sapiens/GRCh38/genome/GRCh38.CHR.fa
W=/media/share/workdir/1610305236/lncProt
mkdir -p $W/data/gbm_own/ref

if [ ! -s $W/data/gbm_own/ref/transcripts.fa ]; then
  singularity exec --no-mount bind-paths --bind $REF:$REF:ro --bind $W:$W \
    /media/share/resource/container/gffread/gffread-0.12.7.sif gffread -w $W/data/gbm_own/ref/transcripts.fa -g $GEN $GTF
fi
echo "transcripts: $(grep -c '^>' $W/data/gbm_own/ref/transcripts.fa)"

if [ ! -d $W/data/gbm_own/ref/salmon_idx ]; then
  singularity exec --no-mount bind-paths --bind $REF:$REF:ro --bind $W:$W \
    /media/share/resource/container/salmon/salmon-1.10.2.sif salmon index \
      -t $W/data/gbm_own/ref/transcripts.fa -i $W/data/gbm_own/ref/salmon_idx \
      -k 31 -p 16
fi
echo "salmon index built: $(date)"
