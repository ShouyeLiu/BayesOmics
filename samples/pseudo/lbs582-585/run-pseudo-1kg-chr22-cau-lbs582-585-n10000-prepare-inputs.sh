#!/bin/sh
# Fixed pseudo10000/cis50k/blocks582-585 website example. No R rerun.
set -eu
example_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_dir=$(CDPATH= cd -- "$example_dir/../../.." && pwd)
binary=${BAYESOMICS_BIN:-$project_dir/build/BayesOmics64}
output=${BAYESOMICS_EXAMPLE_OUTPUT:-$example_dir/output}
threads=${BAYESOMICS_THREADS:-4}
test -x "$binary"
mkdir -p "$output/pseudo-1kg-chr22-cau-lbs582-585-n10000-ld-blocks"
output=$(CDPATH= cd -- "$output" && pwd)
cd "$example_dir/data"
"$binary" --bfile pseudo-1kg-chr22-cau-lbs582-585-n10000 --block-info pseudo-1kg-chr22-cau-lbs582-585-n10000-ref4cM-v37.pos \
  --make-block-ldm --thread "$threads" --out "$output/pseudo-1kg-chr22-cau-lbs582-585-n10000-ld-blocks"
"$binary" --ldm "$output/pseudo-1kg-chr22-cau-lbs582-585-n10000-ld-blocks" --merge-block-ldm-info --out "$output/pseudo-1kg-chr22-cau-lbs582-585-n10000-ld-blocks"
"$binary" --ldm "$output/pseudo-1kg-chr22-cau-lbs582-585-n10000-ld-blocks" --make-ldm-eigen \
  --ldm-eigen-cutoff 1 --thread "$threads" --out "$output/pseudo-1kg-chr22-cau-lbs582-585-n10000-ld-blocks"
"$binary" --bfile pseudo-1kg-chr22-cau-lbs582-585-n10000 --block-info pseudo-1kg-chr22-cau-lbs582-585-n10000-ref4cM-v37.pos \
  --make-eigen --ldm-eigen-cutoff 1 --thread "$threads" --out "$output/pseudo-1kg-chr22-cau-lbs582-585-n10000-ld-direct"
"$binary" --bfile pseudo-1kg-chr22-cau-lbs582-585-n10000 --gene-annotation pseudo-1kg-chr22-cau-lbs582-585-n10000-gene-annotation-chr1-22-hg19-v40.txt \
  --cis-wind 0.05 --make-eigen-gene --ldm-eigen-cutoff 1 --out "$output/pseudo-1kg-chr22-cau-lbs582-585-n10000-annotation-gene-ld"
"$binary" --bfile pseudo-1kg-chr22-cau-lbs582-585-n10000 --beqtl-summary-gz pseudo-1kg-chr22-cau-lbs582-585-n10000-eqtl --make-eigen-gene \
  --ldm-eigen-cutoff 1 --out "$output/pseudo-1kg-chr22-cau-lbs582-585-n10000-molecular-ld"
"$binary" --eqtl-flist pseudo-1kg-chr22-cau-lbs582-585-n10000-eqtl.flist --make-besd --out "$output/pseudo-1kg-chr22-cau-lbs582-585-n10000-eqtl-besd"
"$binary" --beqtl-summary "$output/pseudo-1kg-chr22-cau-lbs582-585-n10000-eqtl-besd" --add-gene-n pseudo-1kg-chr22-cau-lbs582-585-n10000-gene-n.txt \
  --make-query --out "$output/pseudo-1kg-chr22-cau-lbs582-585-n10000-eqtl-roundtrip"
