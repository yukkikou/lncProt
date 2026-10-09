#!/bin/bash
# =============================================================================
# lncProt — lncRNA vs PCG promoter comparative analysis
# Central configuration: ALL paths and parameters (single source of truth).
# Reference tree is READ-ONLY and must never be written to.
# =============================================================================
set -euo pipefail

# ---- Workdir (only writable area) -------------------------------------------
export WORKDIR=/media/share/workdir/1610305236/lncProt
export DATA=$WORKDIR/data
export RESULTS=$WORKDIR/results
export LOGS=$WORKDIR/logs
export FIGS=$WORKDIR/figures

# ---- Read-only reference / database roots ------------------------------------
export REF=/media/share/resource/reference
export DB=/media/share/resource/database
export ANNO=$REF/homo_sapiens/GRCh38/annotation
export GENOME=$REF/homo_sapiens/GRCh38/genome/GRCh38.CHR.fa
export GENCODE=$ANNO/gencode.v46.annotation.gtf
export TE_BED=$ANNO/hg38.nrph.TE.anno.sorted.bed          # Dfam-based TE annotation
export HERV_BED=$ANNO/HERV_noncoding.bed
export CHAIN_HG19_HG38=$REF/homo_sapiens/liftOver/hg19ToHg38.over.chain
export PHASTCONS_BW=$DB/PhyloP_score/hg38.phastCons470way.bw
export PHYLOP_BW=$DB/PhyloP_score/241-mammalian-2020v2.bigWig   # Cactus/Zoonomia 241-mammal phyloP
export REMAP_ALL=$DB/remap2022/remap2022_all_macs2_hg38_v1_0.bed
export REMAP_NR=$DB/remap2022/remap2022_nr_macs2_hg38_v1_0.bed
export TFBS_ATLAS_DIR=$DB/TFBS_GRCh38/PMID37104580          # 368 per-TF BEDs (cross-check)
export HOCOMOCO_LOCAL=$DB/MotifDB/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme
export JASPAR_LOCAL=$DB/JASPAR/JASPAR2024_CORE_non-redundant_pfms_meme.txt

# ---- Compute: python + singularity containers --------------------------------
# NOTE: containers live under /media/share/resource/container (NOT under $REF)
export PY=/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/bin/python
export SIF_BEDTOOLS=/media/share/resource/container/bedtools/bedtools-2.30.0.sif
export SIF_MEME=/media/share/resource/container/meme/meme-5.4.1.sif
export SIF_PYBIGWIG=/media/share/resource/container/pyBigWig/pyBigWig.sif

# ---- Analysis parameters ------------------------------------------------------
export FEAT_UP=500        # promoter window (upstream of TSS) for sequence features
export FEAT_DN=200        # downstream
export FIMO_UP=500        # FIMO window
export FIMO_DN=100
export FIMO_P=1e-4        # FIMO p-value threshold
export DIVERGENT_WIN=1000 # primary divergent-pair cutoff (bp between opposite-strand TSSs)
export DIVERGENT_SENS="500,2000"  # sensitivity analysis cutoffs
export INTERGENIC_MIN=10000       # lncRNA >=10kb from any PCG = far intergenic
export CAGE_MATCH_WIN=250 # CAGE peak summit within +/- bp of annotated TSS
export CAGE_TPM_MIN=1.0   # minimum CAGE TPM for active promoter
export MAXCHUNK_NJOBS=10  # FIMO motif chunks (array tasks)

# ---- Roadmap ChromHMM cell types (hg19 -> liftover) ---------------------------
export CHROMHMM_CELLS="E003 E017 E067 E071 E088"   # H1-ESC, GM12878, liver, brain, testis
export CHROMHMM_URLBASE="https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final"

mkdir -p "$DATA" "$RESULTS" "$RESULTS/classes" "$RESULTS/step_summaries" "$LOGS/failures" "$FIGS" 2>/dev/null || true
