# GenoBridge

[![GitHub release](https://img.shields.io/github/v/release/nps-genomics/GenoBridge)](https://github.com/nps-genomics/GenoBridge/releases)
[![Conda version](https://img.shields.io/conda/vn/nps-genomics/genobridge)](https://anaconda.org/nps-genomics/genobridge)
[![Platform](https://img.shields.io/badge/platform-linux--64-blue)](https://anaconda.org/nps-genomics/genobridge)
[![Python](https://img.shields.io/badge/python-3.12-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-Apache--2.0-green)](LICENSE)

**A sample-size-adaptive machine-learning and mixed-model GWAS framework for genotype-to-phenotype analysis.**

> **Recommended installation:** use the Conda command in the [Installation](#installation) section. The GitHub wheel is provided only as an advanced fallback.

GenoBridge integrates genomic prediction with a predictability-gated genome-wide association study workflow. It:

1. Automatically adjusts model complexity according to sample size.
2. Evaluates multiple prediction models for every phenotype.
3. Identifies traits with sufficient genomic predictability.
4. Runs kinship-adjusted mixed-model GWAS using GEMMA.
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
Genotype VCF ── PLINK2 preprocessing ── GEMMA mixed-model GWAS
                                  │
                                  ├── Association statistics
                                  ├── Manhattan and QQ plots
                                  ├── Genomic inflation estimates
                                  └── Nearest-gene annotation
```

GenoBridge prediction evaluates Ridge regression, Random Forest, XGBoost, neural-network, and stacking models. Model capacity, PCA dimensions, validation strategy, and ensemble composition are adjusted according to dataset size.

The GWAS module uses PLINK2 for genotype preprocessing and GEMMA for kinship-adjusted linear mixed-model association testing.

---
## Software distribution

GenoBridge v1.0.1 is distributed for **Linux x86_64** systems.

The recommended installation method is the public Conda package hosted in the
[`nps-genomics`](https://anaconda.org/nps-genomics/genobridge) channel. This
method creates a clean Python 3.12 environment and installs GenoBridge together
with GEMMA, PLINK2, and all required Python dependencies.

A compiled CPython 3.12 wheel is also provided through GitHub Releases for
advanced users who need a manual installation. The wheel must be installed
inside a Python 3.12 Linux x86_64 environment.

The prediction and GWAS implementations are installed as compiled Linux
extension modules. 

---

## Installation

### Supported systems

GenoBridge v1.0.1 currently supports:

* Linux x86_64 (`linux-64`)
* CPython 3.12
* glibc 2.17 or newer

### Recommended installation

Create a new Conda environment. Do not install GenoBridge into an existing
Python 3.6, 3.10, or 3.11 environment.

```bash
conda create -n genobridge \
  --override-channels \
  -c nps-genomics \
  -c conda-forge \
  -c bioconda \
  genobridge=1.0.1 \
  -y

conda activate genobridge
```

This command installs:

* GenoBridge 1.0.1
* Python 3.12
* GEMMA 0.98.5
* PLINK2
* Required Python libraries

Verify the installation:

```bash
python -c "import genobridge; print('GenoBridge', genobridge.__version__)"

genobridge-predict --help
genobridge-gwas --help
gemma -h
plink2 --version
```

Expected GenoBridge version:

```text
GenoBridge 1.0.1
```

### Advanced manual installation from the GitHub wheel

Most users should use the Conda installation above.

The wheel installation is intended only for users who already understand
Python environments and compiled wheel compatibility. The wheel requires
CPython 3.12 on Linux x86_64.

#### Where to download the wheel

Open the GenoBridge v1.0.1 release page:

<https://github.com/nps-genomics/GenoBridge/releases/tag/v1.0.1>

On the release page, expand **Assets** and download:

```text
genobridge-1.0.1-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.whl
```

Direct wheel download:

<https://github.com/nps-genomics/GenoBridge/releases/download/v1.0.1/genobridge-1.0.1-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.whl>

Do not download GitHub's automatically generated **Source code (zip)** or
**Source code (tar.gz)** files for installation. Those archives contain the
public repository contents, not the compiled GenoBridge package.

Create a clean Python 3.12 environment:

```bash
conda create -n genobridge-wheel \
  --override-channels \
  -c conda-forge \
  -c bioconda \
  python=3.12 \
  pip \
  gemma=0.98.5 \
  plink2 \
  -y

conda activate genobridge-wheel
```

Confirm the Python version before installing:

```bash
python --version
uname -m
```

Expected:

```text
Python 3.12.x
x86_64
```

Install the wheel:

```bash
python -m pip install \
  "https://github.com/nps-genomics/GenoBridge/releases/download/v1.0.1/genobridge-1.0.1-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.whl"
```

Do not run this wheel command from `base` or from an environment using Python
3.6, 3.10, or 3.11.

---

## Quick start

### Step 1: genomic prediction and trait triage

```bash
genobridge-predict \
  --phenotype path/to/phenotype.csv \
  --vcf path/to/genotypes.vcf \
  --accession-col accession_id \
  --output results/ \
  --no-outlier-removal \
  
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

Set the gate to zero to analyze all traits:

```bash
--gate 0
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
* Complete GEMMA association output
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

External programs:

* GEMMA 0.98.5
* PLINK2

Manual wheel users must install GEMMA and PLINK2 separately and ensure that both
programs are available through `PATH`.

---

## Platform availability

The current wheel is:

```text
genobridge-1.0.1-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.whl
```

It is compatible with:

* CPython 3.12
* Linux
* x86_64 processors

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

> Singh, N. P., and Mendu, V. (2026). GenoBridge: an automated sample-size-adaptive machine-learning framework for genomic prediction and predictability-gated genome-wide association analysis. Publication details pending.

The citation should be updated when the manuscript and archived software DOI become available.

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

### GWAS halts with `Only X% of genotype samples matched`

This occurs when the VCF contains many more samples than the phenotype file.
For example, a complete diversity panel may contain hundreds or thousands of
samples, whereas phenotypes may be available for only a subset.

Subset the VCF to the phenotyped accessions before running GenoBridge.

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

The SHA-256 checksum for the v1.0.1 wheel is:

```text
869c9163653cdeb99d09543806dfe92ca21d4e3a9b32dbc17f4d051eb4dd7110
```

Verify a downloaded wheel with:

```bash
sha256sum \
  genobridge-1.0.1-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.whl
```

The calculated checksum must match the value above.

---

## Links

* GitHub repository: <https://github.com/nps-genomics/GenoBridge>
* GitHub releases: <https://github.com/nps-genomics/GenoBridge/releases>
* Conda package: <https://anaconda.org/nps-genomics/genobridge>
