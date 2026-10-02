#!/bin/bash
# The three --phasing-sites variants of one sample at a time, together (8
# threads each, config.yaml:resources), so every variant sees the same machine
# load (the consensus timeouts are wall-clock) and the long consensus tails of
# one variant overlap with the others' work. Then the calls, truvari and the
# Mendelian consistency.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
for s in HG002 HG003 HG004; do
    bash run.sh \
        "$PWD/results/ps_ava/20x/$s/work/svirltile.db" \
        "$PWD/results/ps_tiered/20x/$s/work/svirltile.db" \
        "$PWD/results/ps_ref/20x/$s/work/svirltile.db"
done
bash run.sh stage_benchmark stage_mendel
echo "run_grouped.sh: all done"
