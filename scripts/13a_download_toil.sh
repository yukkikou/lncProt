#!/bin/bash
# =============================================================================
# 13a_download_toil.sh — UCSC Xena Toil matrix (TCGA+GTEx unified RSEM TPM,
# 1.3GB) with self-healing resume loop. Sample typing uses TCGA barcode rules
# and the LOCAL GTEx Sample Attributes file (volume2, read-only).
# =============================================================================
source "$(dirname "$0")/common.sh"
trap 'on_err 13a_toil $BASH_COMMAND' ERR
CURRENT_STEP="13a_toil"
progress_init

TD="$DATA/tissue_disease"
mkdir -p "$TD"
URL="https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz"
TOTAL=1323254426

if [ ! -s "$TD/toil_tpm.gz" ] || [ "$(stat -c%s "$TD/toil_tpm.gz")" != "$TOTAL" ]; then
    echo ">>> downloading Toil TPM matrix ($(date))"
    while true; do
        SZ=$(stat -c%s "$TD/toil_tpm.gz" 2>/dev/null || echo 0)
        [ "$SZ" = "$TOTAL" ] && break
        curl -fsSL --max-time 14400 --speed-limit 1000 --speed-time 120 -C - \
             -o "$TD/toil_tpm.gz" "$URL" || { echo "retry..."; sleep 10; }
    done
    gzip -t "$TD/toil_tpm.gz" && echo "GZIP OK"
fi
echo ">>> toil matrix ready: $(stat -c%s "$TD/toil_tpm.gz") bytes"

# GTEx sample attributes from the read-only volume2 tree (copy, never modify)
SATT="/media/share/data8/volume2/GTEx/phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Attributes.GRU.txt"
if [ -s "$SATT" ]; then
    cp "$SATT" "$TD/gtex_sample_attributes.txt"
    echo ">>> GTEx sample attributes copied ($(wc -l < "$TD/gtex_sample_attributes.txt") lines)"
else
    echo "WARN: GTEx sample attributes missing at $SATT" >&2
fi

# TCGA sample sheet + stage annotations (read-only copy)
cp /media/share/data8/volume2/TCGA/0_data/gdc_sample_sheet.2024-07-09.tsv "$TD/" 2>/dev/null || true
cp /media/share/data8/volume2/TCGA/0_data/stage_anno.txt "$TD/" 2>/dev/null || true
cp /media/share/data8/volume2/TCGA/0_data/clinical.2024.10.20/clinical_simple.tsv "$TD/" 2>/dev/null || true
progress "13a_toil" "DONE" "bytes=$(stat -c%s "$TD/toil_tpm.gz" 2>/dev/null)"
