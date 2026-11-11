#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
stage=${1:-validate}
if [ "$#" -gt 0 ]; then shift; fi
exec python3 "$repo/scripts/run_stage.py" "$stage" "$@"
