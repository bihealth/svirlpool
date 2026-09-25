#!/bin/bash
# usage: run_profile.sh <out_prefix> <profile_batch.py args...>
# Runs profile_batch.py with this checkout's sources and the dev pixi env.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
ENV="${SVIRLPOOL_ENV_PREFIX:-/home/mayv_c/development/svirlpool/.pixi/envs/dev}"
export PATH="${PERF_SHIM:+$PERF_SHIM:}$ENV/bin:$PATH"
export PYTHONPATH="$ROOT/src"
OUT="$1"
shift
mkdir -p "$(dirname "$OUT")"
python "$HERE/profile_batch.py" "${WORKDIR:-/home/mayv_c/development/svp_improvements/results/ava_phase2/20x/HG002/work}" "$OUT" "$@" >"$OUT.out" 2>"$OUT.err"
