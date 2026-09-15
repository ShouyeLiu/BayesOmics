# SBayesCO-EIEO: fixed pseudo example

This example reuses a saved causal simulation. It does not generate new genotypes
or phenotypes, and does not require R to run the C++ analysis.

| Setting | Value |
| --- | --- |
| Individuals | 10,000 pseudo individuals |
| SNPs | 5,092, chromosome 22 |
| LD blocks | 582–585 |
| Simulated molecular phenotypes | 12 |
| Model molecular SNP–gene pairs | 493 |
| Annotation window | Gene midpoint ±50 kb |
| Simulation seed | 20260914 |
| Each model chain | 3,000 iterations; 2,000 burn-in; thinning 1 |
| Model threads | 1, to match the saved reference settings |

The raw annotation map contains 133 genes and 6,539 pairs; the analysis uses the
12 simulated molecules and their 493 pairs. These are distinct selections.

## Run from the repository root

After following the [installation instructions](../../../README.md#installation):

```sh
sh samples/pseudo/lbs582-585/run-pseudo-1kg-chr22-cau-lbs582-585-n10000-prepare-inputs.sh
sh samples/pseudo/lbs582-585/run-pseudo-1kg-chr22-cau-lbs582-585-n10000-SBayesCO-EIEO-individual.sh
sh samples/pseudo/lbs582-585/run-pseudo-1kg-chr22-cau-lbs582-585-n10000-SBayesCO-EIEO-summary.sh
```

The default executable is `build/BayesOmics64`. Set `BAYESOMICS_BIN` to an absolute
path to use another build. `BAYESOMICS_EXAMPLE_OUTPUT` selects an output directory;
its default is this example's `output/`. `BAYESOMICS_THREADS` controls LD preparation
(default 4). Both model scripts deliberately use one thread.

The preparation script builds block LD, block and direct eigen references,
gene LD from annotation and from the supplied molecular summary, and the
flist → BESD → query conversion. The model scripts use the prepared block and
molecular references where required. They resolve molecular file paths from
`data/`, so the supplied `pseudo-1kg-chr22-cau-lbs582-585-n10000-molecular.plist` and `pseudo-1kg-chr22-cau-lbs582-585-n10000-eqtl.flist` remain portable.

## Inputs and outputs

| Input | Purpose |
| --- | --- |
| `data/pseudo-1kg-chr22-cau-lbs582-585-n10000.bed`, `.bim`, `.fam` | Matched pseudo genotypes and identifiers |
| `data/pseudo-1kg-chr22-cau-lbs582-585-n10000-gwas.phen` | Individual GWAS phenotype |
| `data/pseudo-1kg-chr22-cau-lbs582-585-n10000-molecular.plist`, `data/molecular/` | Molecular phenotype manifest and values |
| `data/pseudo-1kg-chr22-cau-lbs582-585-n10000-model-snp-gene-map.txt` | The simulated molecular SNP–gene pairs |
| `data/pseudo-1kg-chr22-cau-lbs582-585-n10000-gwas.ma` | GWAS summary statistics |
| `data/pseudo-1kg-chr22-cau-lbs582-585-n10000-eqtl.query.gz` | Molecular QTL summary statistics |
| `data/pseudo-1kg-chr22-cau-lbs582-585-n10000-genotype-scale.txt` | Empirical genotype scales used for comparison |
| `data/pseudo-1kg-chr22-cau-lbs582-585-n10000-truth.json` | Saved true effects and realized heritabilities |
| `data/pseudo-1kg-chr22-cau-lbs582-585-n10000-input-manifest.json` | Input dimensions, source identifiers and SHA256 hashes |

Both model scripts request complete text MCMC output. Results are written under
`output/models/`, with separate `pseudo-1kg-chr22-cau-lbs582-585-n10000-SBayesCO-EIEO-individual` and `pseudo-1kg-chr22-cau-lbs582-585-n10000-SBayesCO-EIEO-summary` prefixes.
With thinning 1, the posterior summaries retain 1,000 post-burn-in iterations.
Generated LD, logs and MCMC files are local outputs and are excluded from Git.

LD correlation and genotype covariance scales must be aligned. Keep the supplied
`pseudo-1kg-chr22-cau-lbs582-585-n10000-genotype-scale.txt`, allele coding, SNP order and full eigen cutoffs when using
these settings. Per-allele input effects and standardized simulation effects are
not interchangeable without their recorded scales.

For model and file-format explanations, use the
[BayesOmics manual](https://shouyeliu.github.io/softwares/content-softwares.html).
