#!/bin/bash
# Copy of d3/truvari_strict.sh writing into proto/: truvari bench of HG002's
# column vs T2T-Q100 as in bench.sh, plus --pctsize 0.9 --pctseq 0.9.
#   truvari_strict.sh <arm>
set -eu
ARM=$1
ROOT=/home/mayv_c/development/svp_mendel
TRUVARI_BIN=/home/mayv_c/development/svp_merging/.conda-envs/truvari/bin
TRUTH=/home/mayv_c/development/svp_improvements/resources/truthsets
OUT=$ROOT/proto/truvari_strict/$ARM
mkdir -p "$ROOT/proto/truvari_strict"
rm -rf "$OUT"
export PATH="$TRUVARI_BIN:$PATH"
truvari bench -b "$TRUTH/GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz" -c "$ROOT/arms/$ARM/family.vcf.gz" \
    -f "$ROOT/inputs/GRCh38.fa" -o "$OUT" \
    --includebed "$TRUTH/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed" \
    --passonly --pick ac --dup-to-ins --sizemin 50 --cSample HG002 --no-ref c \
    --pctsize 0.9 --pctseq 0.9 > "$OUT.log" 2>&1
python3 -c "import json; s=json.load(open('$OUT/summary.json')); print('$ARM', round(s['precision'],4), round(s['recall'],4), round(s['f1'],4))"
