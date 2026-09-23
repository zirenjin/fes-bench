#!/usr/bin/env bash
set -euo pipefail

# Public files are intentionally excluded from git; hashes and destinations are
# documented in data/README.md and in each data/processed/<system>/download.json.
repo_root="$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)"
destination="$repo_root/data/raw"
mkdir -p "$destination"

download() {
  local url="$1" target="$2" expected="$3"
  mkdir -p "$(dirname "$target")"
  curl --fail --location --continue-at - "$url" --output "$target"
  printf '%s  %s\n' "$expected" "$target" | sha256sum --check --strict
}

if [[ "${1:-}" == "--dry-run" ]]; then
  printf '%s\n' "See data/README.md and data/*/download.json for the configured public URLs and SHA-256 values."
  exit 0
fi

download "https://darus.uni-stuttgart.de/api/access/datafile/380770" "$destination/sio2/direct_upsampling_data_rungs_1-3.zip" "d6767aa491589e6c76d9b929f6d228fc3315e2000fb933f17b554282003bcd79"
download "https://darus.uni-stuttgart.de/api/access/datafile/239818" "$destination/hf/Hf_hcp_PBE.tar.gz" "e79d5c3c6be5bc795181447d654a8cce1f939c25756716bf06cba6be4bbab3b2"
download "https://darus.uni-stuttgart.de/api/access/datafile/239823" "$destination/hf/Hf_bcc_PBE.tar.gz" "8efef64d7d9daf2df5c39718fc09acf9e150212a702d45e5f8d4c61f2e374657"
