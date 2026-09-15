#!/bin/sh
# SBayesCO-EIEO tutorial: fixed pseudo data; same settings as saved 3000/2000 runs.
set -eu
example_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_dir=$(CDPATH= cd -- "$example_dir/../../.." && pwd)
binary=${BAYESOMICS_BIN:-$project_dir/build/BayesOmics64}
output=${BAYESOMICS_EXAMPLE_OUTPUT:-$example_dir/output}
test -x "$binary"
mkdir -p "$output/models"
output=$(CDPATH= cd -- "$output" && pwd)
cd "$example_dir/data"
"$binary" --mcmc-type EIEO --genotype-scale pseudo-1kg-chr22-cau-lbs582-585-n10000-genotype-scale.txt \
  --hsq .5 --hsq-cis .5 \
  --pi-genic .4 --pi-genic-gwas .4 --pi-genic-xqtl .4 --pi-intergenic .01 \
  --chain-length 3000 --burn-in 2000 --thin 1 --seed 20260914 \
  --thread 1 --memory auto --genotype-storage auto \
  --write-mcmc-txt --out-freq 3000 \
  --bayes CO --bfile pseudo-1kg-chr22-cau-lbs582-585-n10000 --pheno pseudo-1kg-chr22-cau-lbs582-585-n10000-gwas.phen \
  --gene pseudo-1kg-chr22-cau-lbs582-585-n10000-molecular.plist --gene-snp-map pseudo-1kg-chr22-cau-lbs582-585-n10000-model-snp-gene-map.txt \
  --out "$output/models/pseudo-1kg-chr22-cau-lbs582-585-n10000-SBayesCO-EIEO-individual"
