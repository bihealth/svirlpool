#!/bin/bash
# Default sv-calling output of the branch (flag off) vs c8804b1 on a chr21/22 crID subset.
set -euo pipefail
ENV=/home/mayv_c/development/svirlpool/.claude/worktrees/strip-error-rate/.pixi/envs/dev
WT=/home/mayv_c/development/svirlpool/.claude/worktrees/agent-a7ab379d77ee1935e
IN=/home/mayv_c/development/svp_mendel/inputs
OUT=/home/mayv_c/development/svp_mendel/proto/default_identity
mkdir -p "$OUT"
ARGS=(--input "$IN/dbs/HG002/svirltile.db" "$IN/dbs/HG003/svirltile.db" "$IN/dbs/HG004/svirltile.db"
      --reference "$IN/GRCh38.fa" --threads 4 --min-sv-size 50 --skip-covtrees
      --candidate-regions-file /home/mayv_c/development/svp_mendel/proto/cr_chr21_22.tsv)
env PATH="$ENV/bin:$PATH" "$ENV/bin/svirlpool" sv-calling "${ARGS[@]}" --output "$OUT/c8804b1.vcf.gz" > "$OUT/c8804b1.log" 2>&1
env PYTHONPATH="$WT/src" PATH="$ENV/bin:$PATH" "$ENV/bin/svirlpool" sv-calling "${ARGS[@]}" --output "$OUT/branch_off.vcf.gz" > "$OUT/branch_off.log" 2>&1
env PYTHONPATH="$WT/src" PATH="$ENV/bin:$PATH" "$ENV/bin/svirlpool" sv-calling "${ARGS[@]}" --output "$OUT/branch_align.vcf.gz" --core-alignment-merge align > "$OUT/branch_align.log" 2>&1
for f in c8804b1 branch_off branch_align; do
    zcat "$OUT/$f.vcf.gz" | grep -v "^##fileDate\|^##source\|^##commandline\|^##cmd" > "$OUT/$f.txt"
done
echo "c8804b1 vs branch_off:"; cmp "$OUT/c8804b1.txt" "$OUT/branch_off.txt" && echo identical
echo "records:"; grep -vc "^#" "$OUT/c8804b1.txt" "$OUT/branch_off.txt" "$OUT/branch_align.txt"
