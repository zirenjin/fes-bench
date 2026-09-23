#!/usr/bin/env bash
set -euo pipefail
repo_root="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
exec bash "$repo_root/src/experiments/data_prep/download.sh" "$@"
