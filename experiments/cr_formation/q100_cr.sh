#!/bin/bash
# HG002 consensuses of the three CR variants vs the Q100 assembly (consensus_q100.py),
# python with edlib from the truvari env (the svirlpool dev env has no edlib),
# minimap2 from the svirlpool dev env.
set -euo pipefail
cd /home/mayv_c/development/svp_tiered15
export PATH=/home/mayv_c/development/svp_tiered15/.conda-envs/truvari/bin:/home/mayv_c/development/svirlpool/.claude/worktrees/strip-error-rate/.pixi/envs/dev/bin:$PATH
Q100=/home/mayv_c/development/svirlpool/.claude/worktrees/ref-vs-ava/experiments/consensus_perf/consensus_q100.py
for v in cr_base cr_s1 cr_s1_sat; do
    python "$Q100" "results/$v/20x/HG002/work" "results/$v/20x/HG002/consensus_q100.tsv" --threads 24
    echo "q100 $v done"
done
