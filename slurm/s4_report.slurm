#!/bin/bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH -J lncP-rep
#SBATCH -p cpu
#SBATCH -D /media/share/workdir/1610305236/lncProt
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -o /media/share/home/1610305236/log/%x-%j.out
#SBATCH -e /media/share/home/1610305236/log/%x-%j.err

source /media/share/home/1610305236/.local/share/mamba/etc/profile.d/mamba.sh
micromamba activate Jupyter

echo ">>> Step 11: figures + report"
python -W ignore scripts/11_report_figures.py || echo "WARN: report incomplete"

source scripts/common.sh
ref_manifest check || echo "FATAL: reference integrity violation"

tar czf "backup/snapshot_s4_$(date +%Y%m%d_%H%M).tar.gz" results/ figures/ REPORT.md 2>/dev/null || true
echo "s4 complete: $(date)"
