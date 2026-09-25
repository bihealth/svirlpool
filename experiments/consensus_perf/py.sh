#!/bin/bash
# usage: py.sh <script.py> [args...]
# Runs a script of this directory with this checkout's sources and the dev pixi env.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
ENV="${SVIRLPOOL_ENV_PREFIX:-/home/mayv_c/development/svirlpool/.pixi/envs/dev}"
export PATH="$ENV/bin:$PATH"
export PYTHONPATH="$ROOT/src"
script="$1"
shift
exec python "$HERE/$script" "$@"
