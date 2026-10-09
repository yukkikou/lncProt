#!/bin/bash
# =============================================================================
# 07b_remap_binding.sh — real TF binding (ReMap 2022 ChIP-seq peaks, hg38).
# The 16GB all-file is STREAMED with bedtools -sorted -g (low memory):
#   1. build genome order file matching the ReMap file's chrom order
#   2. sort promoters with the same order
#   3. stream intersect -> (promoter, TF, cell) presence rows
#   4. nr-file any-TF coverage per promoter (-wo)
# =============================================================================
source "$(dirname "$0")/common.sh"
trap 'on_err 07b_remap $BASH_COMMAND' ERR
CURRENT_STEP="07b_remap"
progress_init

echo ">>> $(date) 07b_remap start"

# ---- genome order file in ReMap chrom order ---------------------------------
GFILE="$DATA/qc/remap_genome.txt"
if [ ! -s "$GFILE" ]; then
    ORDER="$DATA/qc/remap_chrom_order.txt"
    if [ ! -s "$ORDER" ]; then
        awk '!seen[$1]++ {print $1}' "$REMAP_ALL" > "$ORDER"
    fi
    # sizes from the GRCh38 fai; unknown -> large sentinel
    awk 'NR==FNR{sz[$1]=$2; next}{print $1"\t"(($1 in sz)?sz[$1]:250000000)}' \
        <(cut -f1,2 "$REF/homo_sapiens/GRCh38/genome/GRCh38.CHR.fa.fai") "$ORDER" > "$GFILE"
fi
echo "    genome order: $(wc -l < "$GFILE") chroms"

# ---- sort promoters to the same order ----------------------------------------
if [ ! -s "$DATA/promoters_feat.sorted.bed" ]; then
    bt sort -g "$GFILE" -i "$RESULTS/promoters_feat.bed" > "$DATA/promoters_feat.sorted.bed"
fi

OUT1="$DATA/remap_hits.tsv.gz"
if [ ! -s "$OUT1" ]; then
    echo ">>> streaming intersect vs 16GB remap_all ($(date))"
    bt intersect -sorted -g "$GFILE" \
        -a "$DATA/promoters_feat.sorted.bed" -b "$REMAP_ALL" -wa -wb \
      | awk -F'\t' 'BEGIN{OFS="\t"}{n=split($10,a,"."); tf=a[2]; cell=a[3];
            gsub(/_[0-9]+$/,"",cell); print $4,tf,cell,$10}' \
      | sort -u | gzip -c > "$OUT1.part" && mv "$OUT1.part" "$OUT1"
    echo "    remap hits: $(zcat "$OUT1" | wc -l) unique (promoter,TF,cell) rows ($(date))"
else
    echo ">>> remap hits exist, skip"
fi

OUT2="$DATA/remap_nr_cov.tsv"
if [ ! -s "$OUT2" ]; then
    echo ">>> nr-file coverage ($(date))"
    bt intersect -sorted -g "$GFILE" -wo \
        -a "$DATA/promoters_feat.sorted.bed" -b "$REMAP_NR" \
      | awk -F'\t' 'BEGIN{OFS="\t"}{print $4, $13-$12}' \
      > "$DATA/remap_nr_ov.tsv"
    "$PY" - "$DATA/remap_nr_ov.tsv" "$OUT2" <<'PYEOF'
import sys, pandas as pd
df = pd.read_csv(sys.argv[1], sep="\t", names=["gene_id", "ov"])
g = df.groupby("gene_id")["ov"].sum().rename("nr_cov")
g.to_csv(sys.argv[2], sep="\t")
print(f"promoters with any-TF overlap: {len(g)}")
PYEOF
else
    echo ">>> nr coverage exists, skip"
fi

progress "07b_remap" "DONE" "rows=$(zcat "$OUT1" 2>/dev/null | wc -l)"
echo ">>> $(date) 07b_remap complete"
