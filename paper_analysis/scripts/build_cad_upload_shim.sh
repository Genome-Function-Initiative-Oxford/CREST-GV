#!/bin/bash
# Build an in-house-collection-formatted view of the `cad` reference collection
# (symlinks, renamed to the in-house naming convention) so it can be run
# through CREST-GV's in_house_collection_path code path, simulating a user
# re-uploading this collection. Used by fig5_upload_example.py.
set -euo pipefail

DEMO="${1:-$(dirname "$0")/../data/cad_upload_shim}"
SRC=/project/Wellcome_Discovery/datashare/CREST-GV_collection/CAD

mkdir -p "$DEMO/bigwigs" "$DEMO/peaks"

for bw in "$SRC"/bigwigs/*.bw; do
  ln -sf "$bw" "$DEMO/bigwigs/$(basename "$bw")"
done

for pk in "$SRC"/peaks/*_L-tron.bed; do
  name=$(basename "$pk" | sed 's/_L-tron\.bed$//')
  ln -sf "$pk" "$DEMO/peaks/${name}.bed"
done

cp "$SRC/CAD_info.csv" "$DEMO/$(basename "$DEMO")_info.csv"

echo "Built in-house shim at: $DEMO"
