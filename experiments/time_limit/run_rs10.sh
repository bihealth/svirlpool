#!/bin/bash
# The trio over the 10% regions (svp_rs10): both variants of every sample
# (8 threads each, three runs at a time on 24 cores), then calls, truvari,
# Mendelian consistency, then the HG002 consensuses vs the Q100 assembly
# (python with edlib from the truvari env, minimap2 from the svirlpool dev env).
set -euo pipefail
H=/home/mayv_c/development/svp_rs10
WT=/home/mayv_c/development/svirlpool/.claude/worktrees/strip-error-rate
CF=$WT/experiments/time_limit/svp_variants_read_selection.yaml
Q100=/home/mayv_c/development/svirlpool/.claude/worktrees/ref-vs-ava/experiments/consensus_perf/consensus_q100.py
VARIANTS="rs_base rs_sel240"
cd "$H"
targets=()
for s in HG002 HG003 HG004; do
    for v in $VARIANTS; do targets+=("$H/results/$v/20x/$s/work/svirltile.db"); done
done
bash run.sh "${targets[@]}" --configfile "$CF"
echo "run_rs10.sh: runs done $(date -Is)"
bash run.sh stage_benchmark stage_mendel --configfile "$CF"
echo "run_rs10.sh: benchmark and mendel done $(date -Is)"
export PATH=$H/.conda-envs/truvari/bin:$WT/.pixi/envs/dev/bin:$PATH
for v in $VARIANTS; do
    python "$Q100" "results/$v/20x/HG002/work" "results/$v/20x/HG002/consensus_q100.tsv" --threads 24
done
echo "run_rs10.sh: all done $(date -Is)"
