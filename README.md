# GenoBridge

**A sample-size-adaptive machine-learning and mixed-model GWAS framework for genotype-to-phenotype analysis.**

GenoBridge integrates genomic prediction with a predictability-gated genome-wide association study workflow. It:

1. Automatically adjusts model complexity according to sample size.
2. Evaluates multiple prediction models for every phenotype.
3. Identifies traits with sufficient genomic predictability.
4. Runs kinship-adjusted mixed-model GWAS using GEMMA.
5. Annotates significant associations to nearby genes.

## Software distribution

GenoBridge is distributed as a **compiled Python wheel**.

* The original Python and Cython implementation source is not distributed.
* Prediction and GWAS modules are installed as compiled Linux extension modules.
* The current release supports:

  * Linux x86_64
  * CPython 3.12
* GEMMA and PLINK2 are required external programs.

The public GitHub repository contains documentation, example input formats, citation information, and release files. It does not contain the private implementation source.

---

## Workflow

```text
Phenotype CSV ┐
              ├── GenoBridge prediction
Genotype VCF ┘         │
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

## Installation

### System compatibility

The current GenoBridge wheel requires:

* Linux x86_64
* Python 3.12
* GEMMA 0.98.5 or later
* PLINK2

### Recommended installation

Create a Conda environment containing Python, GEMMA, and PLINK2:

```bash
conda create -y -n genobridge \
  -c conda-forge \
  -c bioconda \
  python=3.12 \
  pip \
  gemma=0.98.5 \
  plink2

conda activate genobridge
```

Install the compiled GenoBridge wheel:

```bash
python -m pip install   https://github.com/nps-genomics/GenoBridge/releases/download/v1.0.0/ \ genobridge-1.0.0-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.whl

```

Alternatively, download the wheel from the GitHub Releases page and install it locally:

```bash
python -m pip install \
  genobridge-1.0.0-cp312-cp312-linux_x86_64.whl
```

Verify the installation:

```bash
genobridge-predict --help
genobridge-gwas --help
plink2 --version
gemma -h
```

---

## Quick start

### Step 1: genomic prediction and trait triage

```bash
genobridge-predict \
  --phenotype examples/rice_phenotype.csv \
  --vcf examples/rice_geno.vcf \
  --accession-col accession_id \
  --output results/ \
  --no-outlier-removal \
  --max-missing-frac 0.31
```

The rice benchmark uses `0.31` to retain a phenotype with slightly more than 30% missing observations. For general analyses, the default `0.30` threshold is recommended.

### Step 2: mixed-model GWAS

```bash
genobridge-gwas \
  --label rice \
  --phenotype examples/rice_phenotype.csv \
  --vcf examples/rice_geno.vcf \
  --gff examples/rice_genes.gff3 \
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

Python dependencies are installed automatically with the wheel:

* NumPy
* pandas
* SciPy
* scikit-learn
* Matplotlib
* XGBoost
* scikit-allel

External programs must be installed separately:

* GEMMA
* PLINK2

---

## Platform availability

The current release file:

```text
genobridge-1.0.0-cp312-cp312-linux_x86_64.whl
```

is compatible only with:

* CPython 3.12
* Linux
* x86_64 processors

Additional wheels are required for other Python versions, operating systems, or CPU architectures.

---

## License

GenoBridge is proprietary research software.

Copyright © 2026 Nagendra Pratap Singh. All rights reserved.

The software is distributed in compiled form for authorized use. The original implementation source is not distributed. Redistribution, modification, sublicensing, reverse engineering, and commercial use are prohibited unless expressly authorized in writing by the copyright holder.

The final license should be reviewed and approved according to applicable institutional intellectual-property requirements.

---

## Citation

Users of GenoBridge in academic research should cite:

> Singh, N. P., and Mendu, V. (2026). GenoBridge: an automated sample-size-adaptive machine-learning framework for genomic prediction and predictability-gated genome-wide association analysis. Publication details pending.

A machine-readable citation is provided in:

```text
CITATION.cff
```

The archived software release DOI will be added after release deposition.

---

## Troubleshooting

### `gemma: command not found`

Install GEMMA and make sure it is available through `PATH`:

```bash
conda install -c conda-forge -c bioconda gemma=0.98.5
```

### `plink2: command not found`

Install PLINK2:

```bash
conda install -c conda-forge -c bioconda plink2
```

### `--max-missing-frac: expected one argument`

The option requires a numeric value:

```bash
--max-missing-frac 0.30
```

### Wheel is not supported on this platform

Confirm that the environment uses CPython 3.12 on Linux x86_64:

```bash
python --version
uname -m
```

### Large VCF consumes too much memory

LD-prune the VCF before prediction:

```bash
plink2 \
  --vcf input.vcf \
  --indep-pairwise 100 10 0.2 \
  --out pruning
```

The full VCF may still be used in the GWAS stage.

---

## Release integrity

A SHA-256 checksum is provided for each wheel release.

Verify a downloaded wheel with:

```bash
sha256sum genobridge-1.0.0-cp312-cp312-linux_x86_64.whl
```
