#!/usr/bin/env bash
# Import a list of files verbatim from an upstream checkout into extern/<name>.
# Usage: import-files.sh <git-checkout> <name> <file-list>
# The file list holds one path per line, relative to the checkout.
set -euo pipefail
src=$1; name=$2; list=$3
root=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
dst=$root/extern/$name
rm -rf "$dst"; mkdir -p "$dst"
git -C "$src" archive --format=tar HEAD -- $(sort -u "$list") | tar -x -C "$dst"
