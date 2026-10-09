#!/bin/bash
# =============================================================================
# run_all.sh — submit the full DAG with dependency chaining.
#   s0 download -> s1 prepare -> {s2a..s2e} parallel -> s3 enrichment -> s4 report
# All stages restartable (skip-if-exists), so rerunning is safe.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."

S0=$(sbatch --parsable slurm/s0_download.slurm)
echo "s0 download  -> $S0"
S1=$(sbatch --parsable --depend=afterok:$S0 slurm/s1_prepare.slurm)
echo "s1 prepare   -> $S1"
S2A=$(sbatch --parsable --depend=afterok:$S1 slurm/s2a_seq_te.slurm)
S2B=$(sbatch --parsable --depend=afterok:$S1 slurm/s2b_cage.slurm)
S2C=$(sbatch --parsable --depend=afterok:$S1 slurm/s2c_consv.slurm)
S2D=$(sbatch --parsable --depend=afterok:$S1 slurm/s2d_fimo.slurm)
S2E=$(sbatch --parsable --depend=afterok:$S1 slurm/s2e_remap.slurm)
echo "s2 parallel  -> $S2A $S2B $S2C $S2D(array) $S2E"
S3=$(sbatch --parsable --depend=afterok:$S2A:$S2B:$S2C:$S2D:$S2E slurm/s3_enrichment.slurm)
echo "s3 enrich    -> $S3"
S4=$(sbatch --parsable --depend=afterok:$S3 slurm/s4_report.slurm)
echo "s4 report    -> $S4"
echo "DAG submitted. Monitor: squeue -u \$USER ; logs in ~/log/"
