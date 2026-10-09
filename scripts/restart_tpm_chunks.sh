#!/bin/bash
# restart_tpm_chunks.sh — resume stalled parallel range downloads of the
# FANTOM TPM matrix with append-mode loops (self-healing, stall-killing).
cd /media/share/workdir/1610305236/lncProt
TOTAL=844538065
NCHUNK=6
CSIZE=$(( (TOTAL + NCHUNK - 1) / NCHUNK ))
URL="https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz"
for i in 0 1 2 3 5; do
  BASE=$((i * CSIZE))
  END=$(( (i + 1) * CSIZE - 1 ))
  [ $END -ge $TOTAL ] && END=$((TOTAL - 1))
  WANT=$((END - BASE + 1))
  CUR=$(stat -c%s data/fantom/par/chunk_$i 2>/dev/null || echo 0)
  if [ "$CUR" -ge "$WANT" ]; then
    echo "chunk_$i complete ($CUR)"
    continue
  fi
  echo "restarting chunk_$i at $CUR/$WANT"
  setsid nohup bash -c "
    while true; do
      SZ=\$(stat -c%s data/fantom/par/chunk_$i 2>/dev/null || echo 0)
      [ \$SZ -ge $WANT ] && break
      curl -fsSL --max-time 7200 --speed-limit 500 --speed-time 60 \
           -r \$(( $BASE + \$SZ ))-$END '$URL' >> data/fantom/par/chunk_$i
      sleep 5
    done" >> logs/par_restart_$i.log 2>&1 &
done
sleep 2
pgrep -c curl || true
