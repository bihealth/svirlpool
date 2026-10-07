#!/bin/bash
# Continuation of run_grouped.sh (stopped while HG002 was in its consensus
# tail, its snakemake left running): the three variants of the next sample
# start together once every variant of the previous sample has <= TAIL
# consensus batches left, so the long tails no longer leave cores idle. Then
# the calls, truvari and the Mendelian consistency once nothing else runs.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
TAIL=3
VARIANTS="ps_ava ps_tiered ps_ref"

in_tail() {  # every variant of sample $1 has <= TAIL batches left (or its DB)
    local s=$1 v w n done_
    for v in $VARIANTS; do
        w="$PWD/results/$v/20x/$s/work"
        [[ -s "$w/svirltile.db" ]] && continue
        [[ -s "$w/consensus_batches.tsv" ]] || return 1
        n=$(($(wc -l < "$w/consensus_batches.tsv") - 1))
        done_=$(ls "$w/benchmarks/consensus" 2> /dev/null | wc -l)
        ((n - done_ <= TAIL)) || return 1
    done
}

start_sample() {
    local s=$1
    echo "run_rest.sh: $(date -Is) starting $s"
    bash run.sh \
        "$PWD/results/ps_ava/20x/$s/work/svirltile.db" \
        "$PWD/results/ps_tiered/20x/$s/work/svirltile.db" \
        "$PWD/results/ps_ref/20x/$s/work/svirltile.db" > "logs_run_$s.txt" 2>&1 &
    PIDS+=($!)
}

PIDS=()
prev=HG002
for s in HG003 HG004; do
    until in_tail "$prev"; do sleep 60; done
    start_sample "$s"
    prev=$s
done
# the HG002 snakemake of run_grouped.sh, then the sample runs started here
while pgrep -f "svp_tiered15/results/ps_ava/20x/HG002/work/svirltile.db" > /dev/null; do sleep 60; done
fail=0
for p in "${PIDS[@]}"; do wait "$p" || fail=1; done
echo "run_rest.sh: $(date -Is) runs done (fail=$fail)"
bash run.sh stage_benchmark stage_mendel
echo "run_rest.sh: all done"
