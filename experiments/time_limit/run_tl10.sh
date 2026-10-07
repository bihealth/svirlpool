#!/bin/bash
# HG002 over the 10% regions (svp_tl10): the four time-limit variants together
# (6 threads each, 24 cores), then calls and truvari, then the consensuses vs
# the Q100 assembly (python with edlib from the truvari env, minimap2 from the
# svirlpool dev env).
set -euo pipefail
H=/home/mayv_c/development/svp_tl10
WT=/home/mayv_c/development/svirlpool/.claude/worktrees/strip-error-rate
CF=$WT/experiments/time_limit/svp_variants_time_limit.yaml
Q100=/home/mayv_c/development/svirlpool/.claude/worktrees/ref-vs-ava/experiments/consensus_perf/consensus_q100.py
VARIANTS="tl_off tl_120 tl_180 tl_300"
cd "$H"
targets=()
for v in $VARIANTS; do targets+=("$H/results/$v/20x/HG002/work/svirltile.db"); done
bash run.sh "${targets[@]}" --configfile "$CF"
echo "run_tl10.sh: runs done"
bash run.sh stage_benchmark --configfile "$CF"
echo "run_tl10.sh: benchmark done"
export PATH=$H/.conda-envs/truvari/bin:$WT/.pixi/envs/dev/bin:$PATH
for v in $VARIANTS; do
    python "$Q100" "results/$v/20x/HG002/work" "results/$v/20x/HG002/consensus_q100.tsv" --threads 24
done
echo "run_tl10.sh: all done"
