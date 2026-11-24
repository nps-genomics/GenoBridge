# GenoBridge

[![GitHub release](https://img.shields.io/github/v/release/nps-genomics/GenoBridge)](https://github.com/nps-genomics/GenoBridge/releases)
[![PyPI version](https://img.shields.io/pypi/v/genobridge)](https://pypi.org/project/genobridge/)
[![Conda version](https://img.shields.io/conda/vn/nps-genomics/genobridge)](https://anaconda.org/nps-genomics/genobridge)
[![Platform](https://img.shields.io/badge/platform-linux--64-blue)](https://anaconda.org/nps-genomics/genobridge)
[![Python](https://img.shields.io/badge/python-3.12-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-Apache--2.0-green)](LICENSE)

**GenoBridge: A sample-size-adaptive machine-learning pipeline for genomic prediction and predictability-gated mixed-model GWAS.**

> **Quick installation:** `python -m pip install genobridge` (CPython 3.12, Linux x86-64). For GWAS, install **GEMMA** and **PLINK2** separately, or use the Conda/Micromamba distribution that provides them.

GenoBridge integrates genomic prediction with a predictability-gated genome-wide association study workflow. It:

1. Automatically adjusts model complexity according to sample size.
2. Evaluates multiple prediction models for every phenotype.
3. Identifies traits with sufficient genomic predictability.
4. Runs kinship-adjusted mixed-model GWAS.
5. Annotates significant associations to nearby genes.


---

## Workflow

```text
Phenotype CSV ┐
              ├── GenoBridge prediction
Genotype VCF  ┘        │
                       ├── Per-trait prediction metrics
                       ├── Best model selection
                       └── Predictability gate
                                  │
                                  ▼
Genotype VCF ──  preprocessing ──  mixed-model GWAS
                                  │
                                  ├── Association statistics
                                  ├── Manhattan and QQ plots
                                  ├── Genomic inflation estimates
                                  └── Nearest-gene annotation
```

GenoBridge prediction evaluates Ridge regression, Random Forest, XGBoost, neural-network, and stacking models. Model capacity, PCA dimensions, validation strategy, and ensemble composition are adjusted according to dataset size.


---
## Software distribution

GenoBridge 1.0.1 is available through both **PyPI** and the **`nps-genomics` Conda channel**.
The prediction and GWAS implementations are distributed as compiled extension modules;
the original analysis `.py` and Cython `.pyx` files are **not** included in the published binary wheel.
Compiled binaries are not immune to reverse engineering.

- [PyPI package](https://pypi.org/project/genobridge/1.0.1/)
- [Conda package](https://anaconda.org/nps-genomics/genobridge)
- [GitHub Releases](https://github.com/nps-genomics/GenoBridge/releases/tag/v1.0.1)

## Installation

### Supported platform

The currently published PyPI binary wheel targets **CPython 3.12, Linux x86-64,
glibc 2.17 or newer**. The wheel is not built for macOS, Windows, or other Python versions.
The Python dependencies selected by pip may have stricter platform requirements.

### Option 1: pip (simple prediction installation)

Create a clean Python 3.12 environment, then install the published package:

```bash
python3.12 -m venv genobridge-env
source genobridge-env/bin/activate
python -m pip install --upgrade pip
python -m pip install genobridge==1.0.1
```

Verify the installed commands:

```bash
python -m pip check
genobridge-predict --help
genobridge-gwas --help
```

**Important:** the PyPI wheel installs the GenoBridge Python package and its
declared Python dependencies, **not** the external GEMMA and PLINK2 executables.
The GWAS command requires both on `PATH`. Installing with pip alone is not
sufficient to execute a full GWAS.

For example, in a separate Conda environment with appropriate channels:

```bash
conda create -n genobridge-pip -c conda-forge -c bioconda \
  python=3.12 pip gemma plink2 -y
conda activate genobridge-pip
python -m pip install genobridge==1.0.1
command -v gemma
command -v plink2
```

Use the same environment for pip, GEMMA, and PLINK2 so the executables can be
discovered during GWAS. Check available GEMMA versions against the requirement
of **0.98.5 or newer**.

### Option 2: Micromamba (recommended for complete GWAS setup)

This approach is useful on Linux servers and HPC systems and uses the existing
Conda distribution.

```bash
mkdir -p "$HOME/micromamba-bin"
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
  | tar -xvj -C "$HOME/micromamba-bin" bin/micromamba

export MAMBA_ROOT_PREFIX="$HOME/micromamba"
"$HOME/micromamba-bin/bin/micromamba" create \
  -n genobridge \
  --override-channels \
  -c nps-genomics -c conda-forge -c bioconda \
  python=3.12 genobridge=1.0.1 -y

"$HOME/micromamba-bin/bin/micromamba" run \
  -n genobridge genobridge-predict --help
"$HOME/micromamba-bin/bin/micromamba" run \
  -n genobridge genobridge-gwas --help
```

### Option 3: Conda

```bash
conda create -n genobridge \
  --override-channels \
  -c nps-genomics -c conda-forge -c bioconda \
  python=3.12 genobridge=1.0.1 -y

conda activate genobridge
genobridge-predict --help
genobridge-gwas --help
```

For older Conda installations with slow dependency resolution, Micromamba
may be more convenient.

### Option 4: GitHub wheel (advanced)

The GitHub [v1.0.1 release](https://github.com/nps-genomics/GenoBridge/releases/tag/v1.0.1)
also provides a prebuilt wheel. For most users, **pip from PyPI** is simpler
than downloading and installing a wheel manually.

## Software Registries

GenoBridge is registered in the ELIXIR bioinformatics
software registry, bio.tools.

- bio.tools ID: biotools:genobridge
- https://bio.tools/genobridge

### Installation troubleshooting

- `genobridge-predict: command not found`: activate the environment where
  GenoBridge was installed, then run `python -m pip show genobridge`.
- `gemma: command not found` or `plink2: command not found`: install the
  external executable and make sure it is on `PATH` in the active environment.
- `No matching distribution found`: verify CPython 3.12, Linux x86-64, and
  compatible glibc and Python dependency wheels.
- `FileNotFoundError` for `phenotype.csv` or `genotype.vcf`: replace these
  example filenames with actual input file paths.
- CLI `--help` success checks installation only; full scientific validation
  requires running real data and comparing validated outputs.

---

## Quick start

### Step 1: genomic prediction and trait triage

```bash
genobridge-predict \
  --phenotype path/to/phenotype.csv \
  --vcf path/to/genotypes.vcf \
  --accession-col accession_id \
  --output results/ \
  --no-outlier-removal
```

For general analyses, the default `--max-missing-frac 0.30` threshold is recommended.

### Step 2: mixed-model GWAS

```bash
genobridge-gwas \
  --label rice \
  --phenotype path/to/phenotype.csv \
  --vcf path/to/genotypes.vcf \
  --gff path/to/genes.gff3 \
  --ml-results results/phenotype_prediction_results.csv \
  --gate 0.30 \
  --output results/gwas/
```

---

## Tool 1: genomic prediction

```bash
genobridge-predict \
  --phenotype FILE \
  --vcf FILE \
  [options]
```

| Argument               |        Default | Description                                  |
| ---------------------- | -------------: | -------------------------------------------- |
| `--phenotype`          |       required | Phenotype CSV                                |
| `--vcf`                |       required | Genotype VCF                                 |
| `--accession-col`      | `accession_id` | Accession identifier column                  |
| `--exclude-cols`       | `replicate_id` | Non-phenotype columns to exclude             |
| `--output`             |     `results/` | Output directory                             |
| `--no-outlier-removal` |            off | Disable per-trait outlier removal            |
| `--max-missing-frac`   |         `0.30` | Maximum permitted phenotype missing fraction |
| `--n-components`       |      automatic | Override the number of PCA components        |
| `--test-size`          |      automatic | Override the test-set fraction               |
| `--force-regime`       |          unset | Advanced sample-size regime override         |

Traits exceeding `--max-missing-frac` are excluded and recorded in:

```text
dropped_traits_report.csv
```

The value must be supplied as a fraction:

```bash
--max-missing-frac 0.30
```

Valid examples include:

```text
0.10 = 10%
0.30 = 30%
0.31 = 31%
0.50 = 50%
```

Using `--max-missing-frac 1` permits traits with up to 100% missingness and is not recommended for normal analysis.

---

## Tool 2: mixed-model GWAS

```bash
genobridge-gwas \
  --phenotype FILE \
  --vcf FILE \
  --gff FILE \
  [options]
```

| Argument             |         Default | Description                                      |
| -------------------- | --------------: | ------------------------------------------------ |
| `--phenotype`        |        required | Phenotype CSV                                    |
| `--vcf`              |        required | Genotype VCF                                     |
| `--gff`              |        required | Gene annotation in GFF3 format                   |
| `--ml-results`       |            none | Prediction results used for trait gating         |
| `--gate`             |          `0.30` | Minimum prediction correlation required for GWAS |
| `--gene-function`    |            none | Optional `gene_id,function` CSV                  |
| `--accession-col`    |  `accession_id` | Accession identifier column                      |
| `--exclude-cols`     |  `replicate_id` | Non-phenotype columns                            |
| `--maf`              |          `0.05` | Minor-allele-frequency threshold                 |
| `--geno-missing`     |          `0.10` | Per-marker missingness threshold                 |
| `--bonferroni-alpha` | program default | Bonferroni significance level                    |
| `--suggestive`       | program default | Suggestive threshold setting                     |
| `--label`            |        optional | Cosmetic run label                               |
| `--output`           | program default | GWAS output directory                            |

The `--gate` threshold is applied to the prediction metric used by the
GWAS trait-selection procedure. Its default is `0.30`. For a permissive gate,
check the behavior against the reported prediction metrics and use `--gate -1`
only when it matches the intended inclusion rule.

```bash
genobridge-gwas --phenotype phenotype.csv --vcf genotype.vcf \
  --gff annotation.gff3 --gate 0.30 --output results/gwas/
```

---

## Input formats

### Phenotype CSV

The phenotype file must contain one accession identifier column and one or more phenotype columns.

```csv
accession_id,Brix,Amylose,Plant_height
S001,18.2,22.1,95
S002,17.5,,102
S003,,20.8,88
```

Missing values may be blank or represented as `NA`.

### Genotype VCF

A standard plain-text or compressed VCF may be used. VCF sample identifiers must correspond to phenotype accession identifiers.

Only biallelic SNPs are used for prediction and GWAS preprocessing.

### GFF3 annotation

The GWAS module requires a gene annotation file for nearest-gene assignment.

Functional descriptions may be extracted from attributes such as:

```text
Note=
description=
product=
```

### Optional gene-function table

An optional CSV can provide or override gene descriptions:

```csv
gene_id,function
Gene001,Transcription factor
Gene002,Disease resistance protein
```

---

## Prediction outputs

The prediction output directory contains files including:

### `phenotype_prediction_results.csv`

Per-trait model evaluation results, including:

* Pearson correlation
* Spearman correlation
* Root mean squared error
* Cross-validated R²
* Best-performing model
* Best prediction correlation

This file is supplied to `genobridge-gwas` using `--ml-results`.

### `predictions.csv`

Observed and predicted phenotype values for the held-out test samples.

### `dropped_traits_report.csv`

Traits excluded because their missing-data fraction exceeded `--max-missing-frac`.

### Figures

Prediction figures include model-comparison plots, prediction correlations, and observed-versus-predicted plots.

---

## GWAS outputs

The GWAS output directory contains:

* `gwas_summary.csv`
* Significant and suggestive SNP tables
* Manhattan plots
* QQ plots
* Per-trait genomic inflation estimates
* Nearest-gene annotation
* Optional gene-function annotation

A significant-SNP table may contain:

```text
Chromosome
SNP_ID
Position_bp
Effect_allele
Ref_allele
Allele_freq
Effect_size
SE
P_value
Neg_log10_P
Lambda_REML
Nearest_gene
Distance_bp
Gene_function
```

---

## Runtime dependencies

The recommended Conda installation resolves all declared dependencies
automatically.

Principal Python dependencies include:

* NumPy
* pandas
* SciPy
* scikit-learn
* Matplotlib
* XGBoost
* scikit-allel

Manual wheel users must install GEMMA and PLINK2 separately and ensure that both
programs are available through `PATH`.

---

## Platform availability

The current wheel is:

```text
genobridge-1.0.1-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.whl
```

It is compatible with OS:


* Linux
* Centos
* Ubuntu

The current Conda package is:

```text
genobridge-1.0.1-py312_0
```

and is available for:

```text
linux-64
```

Additional builds are required for other Python versions, operating systems, and CPU architectures.

---

## License

GenoBridge is licensed under the **Apache License 2.0**.

The distributed prediction and GWAS modules are compiled extension modules, and the original implementation source is not included in this public repository. Use and redistribution remain subject to the terms in [`LICENSE`](LICENSE).

Copyright © 2026.

---

## Citation

Users of GenoBridge in academic research should cite:

> Singh, N.P., Mendu, V. GenoBridge: a sample-size-adaptive machine-learning pipeline for genomic prediction and predictability-gated mixed-model GWAS. Plant Methods (2026). https://doi.org/10.1186/s13007-026-01596-5.


> Competing interests:
The authors declare that a patent application related to the technology/method described in this manuscript is pending.

---

## Troubleshooting

### `PackagesNotFoundError` or HTTP 404 during Conda installation

First confirm that the public channel is reachable:

```bash
curl -I -L \
  https://conda.anaconda.org/nps-genomics/linux-64/repodata.json
```

Then use the recommended installation command exactly as shown:

```bash
conda create -n genobridge \
  --override-channels \
  -c nps-genomics \
  -c conda-forge \
  -c bioconda \
  genobridge=1.0.1 \
  -y
```

An HTTP 404 commonly indicates an outdated Conda client, institutional mirror,
proxy, malformed channel URL, or copied non-standard spaces. It does not mean
that the public package is private.

### GenoBridge is incompatible with the existing Python installation

GenoBridge v1.0.1 is compiled for CPython 3.12. Installation into an
environment using Python 3.6, 3.10, or 3.11 will fail with an
`UnsatisfiableError`.

Check the active Python version:

```bash
python --version
```

Create a separate GenoBridge environment instead of upgrading the Python
version of an existing analysis environment.

### Wheel is not supported on this platform

The wheel requires:

* CPython 3.12
* Linux
* x86_64 architecture

Check the active environment:

```bash
command -v python
python --version
uname -m
```

The Python version must be `3.12.x`, and the architecture must be `x86_64`.

### `gemma: command not found`

Activate the recommended environment:

```bash
conda activate genobridge
command -v gemma
```

Manual wheel users can install GEMMA with:

```bash
conda install -c conda-forge -c bioconda gemma=0.98.5
```

### `plink2: command not found`

Activate the recommended environment:

```bash
conda activate genobridge
command -v plink2
```

Manual wheel users can install PLINK2 with:

```bash
conda install -c conda-forge -c bioconda plink2
```

### `--max-missing-frac: expected one argument`

The option requires a numeric fraction:

```bash
--max-missing-frac 0.30
```

### GWAS halts with `Only X% of genotype samples matched` or analyzes too few SNPs

Two different problems can produce an incorrect or interrupted GWAS run.

#### 1. The LD-pruned prediction VCF was supplied to GWAS

The prediction and GWAS stages may use different marker sets. An LD-pruned VCF
can be used for genomic prediction to reduce dimensionality and computational
cost, but GWAS should normally use the full quality-controlled marker set.

Do not run GWAS on the pruned prediction VCF unless that is the intended
analysis. Using a heavily pruned file can substantially reduce genomic
coverage and association power.

For example:

```bash
# Prediction: LD-pruned markers
genobridge-predict \
  --phenotype phenotype.csv \
  --vcf genotype_pruned.vcf \
  --output results/prediction/

# GWAS: full quality-controlled markers
genobridge-gwas \
  --phenotype phenotype.csv \
  --vcf genotype_full_matched.vcf.gz \
  --gff annotation.gff3 \
  --ml-results results/prediction/phenotype_prediction_results.csv \
  --gate 0.30 \
  --output results/gwas/
```

Check the GenoBridge and PLINK2 logs to confirm that the number of analyzed
variants is consistent with the full GWAS marker set rather than the smaller
pruned prediction set.

#### 2. All phenotyped accessions matched, but the VCF contains many additional samples

GenoBridge may report a low overall match percentage when the VCF contains a
large reference panel but phenotypes are available for only a subset. For
example, a message such as:

```text
10,955/44,429 samples matched (24.7%)
```

can still mean that all 10,955 phenotyped accessions were found. The percentage
is low because it is calculated relative to the much larger VCF sample set.

The recommended solution is to subset the **full quality-controlled VCF** to
the phenotyped accessions before running GWAS. Do not subset the already
LD-pruned prediction VCF.

```bash
FULL_VCF=genotype_full.vcf
PHENO=phenotype.csv
PREFIX=genotype_full_matched

# Extract VCF sample IDs.
if [[ "$FULL_VCF" == *.gz ]]; then
  zgrep -m1 '^#CHROM' "$FULL_VCF"
else
  grep -m1 '^#CHROM' "$FULL_VCF"
fi | cut -f10- | tr '\t' '\n' | tr -d '\r' | sort -u > vcf_samples.txt

# Extract phenotype accession IDs from the first CSV column.
awk -F',' 'NR > 1 {gsub(/\r/, "", $1); print $1}' "$PHENO" \
  | sort -u > phenotype_ids.txt

# Keep IDs present in both files.
grep -Fxf phenotype_ids.txt vcf_samples.txt > matched_samples.iid

echo "Phenotype IDs:"
wc -l phenotype_ids.txt

echo "Matched IDs:"
wc -l matched_samples.iid

# Subset the full VCF by sample while retaining its full marker set.
plink2 \
  --vcf "$FULL_VCF" \
  --keep matched_samples.iid \
  --export vcf bgz \
  --out "$PREFIX" \
  --allow-extra-chr
```

PLINK2 accepts a one-column `--keep` file as a list of individual IDs. The
resulting file is:

```text
genotype_full_matched.vcf.gz
```

Confirm the retained sample and marker counts:

```bash
# Number of samples
zgrep -m1 '^#CHROM' genotype_full_matched.vcf.gz \
  | awk -F'\t' '{print NF - 9}'

# Number of VCF variant records
zgrep -vc '^#' genotype_full_matched.vcf.gz
```

Then run GenoBridge GWAS using the sample-matched **full-marker** VCF:

```bash
genobridge-gwas \
  --label study_name \
  --phenotype phenotype.csv \
  --vcf genotype_full_matched.vcf.gz \
  --gff annotation.gff3 \
  --ml-results results/prediction/phenotype_prediction_results.csv \
  --output results/gwas/
```

For large cohorts, GEMMA kinship-matrix construction and eigendecomposition can
require substantial memory and runtime. Also confirm that chromosome names in
the VCF and GFF3 files use compatible conventions, such as `1A` versus
`chr1A`, if nearest-gene annotations are missing.

### Large VCF consumes too much memory

LD-prune the VCF before prediction:

```bash
plink2 \
  --vcf input.vcf \
  --indep-pairwise 100 10 0.2 \
  --out pruning
```

The full VCF may still be used during the GWAS stage.

---

## Release integrity

When distributing a wheel downloaded from GitHub or PyPI, verify its
SHA-256 checksum against the digest published for that **exact file**.
The GitHub release asset and PyPI wheel may not have identical bytes, so
do not reuse a checksum across distributions without verifying it.

```bash
sha256sum genobridge-1.0.1-*.whl
```

---

## Links

* GitHub repository: <https://github.com/nps-genomics/GenoBridge>
* GitHub releases: <https://github.com/nps-genomics/GenoBridge/releases>
* PyPI: <https://pypi.org/project/genobridge/>
* Conda package: <https://anaconda.org/nps-genomics/genobridge>
