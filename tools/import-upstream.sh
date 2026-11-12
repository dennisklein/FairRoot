#!/usr/bin/env bash
# Import one upstream source tree verbatim into extern/<name>.
# Usage: import-upstream.sh <git-checkout> <name> [excluded-top-level-path...]
set -euo pipefail
src=$1; name=$2; shift 2
root=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
dst=$root/extern/$name
spec=(.)
for e in "$@"; do spec+=(":(exclude)$e"); done
rm -rf "$dst"; mkdir -p "$dst"
git -C "$src" archive --format=tar HEAD -- "${spec[@]}" | tar -x -C "$dst"
