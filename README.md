# GenoBridge

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20CentOS%2FRHEL-blue)
![Version](https://img.shields.io/badge/version-1.0-green)

**Adaptive ensemble ML pipeline for genotype-to-phenotype prediction with integrated GWAS and candidate gene discovery.**

Genome-wide association studies (GWAS) utilize genomic information to identify candidate loci associated with agronomically important traits, while genomic prediction (GP) estimates the overall genetic contribution to phenotypic variation. Together, GP and GWAS are foundational methodologies in modern plant breeding for crop improvement. However, none of the existing pipelines integrates all three components, sample-size-adaptive ML model selection, CV-R² evaluation, and ML-gated GWAS within a single automated workflow. Current analytical methods exhibit three major limitations, model selection is largely assumption-driven, overlooking the predictable relationship between cohort size and optimal model class, Pearson correlation remains the default accuracy metric despite its inability to distinguish true predictive generalization from rank-order preservation; and GWAS and GP are traditionally conducted independently, misallocating statistical power by running association analysis on traits with no evidence of genomic predictability.  



Supports 30–100,000+ samples across any plant species. No Python installation required — runs as a standalone binary.

---

## Table of Contents

- [Overview](#overview)
- [Tools Included](#tools-included)
- [Download](#download)
- [Installation](#installation)
- [Usage](#usage)
- [Input File Formats](#input-file-formats)
- [Output Files](#output-files)
- [System Requirements](#system-requirements)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Overview

GenoBridge is a two-tool pipeline designed for plant genomics research:

1. **GenoBridge** — Adaptive ensemble machine learning for genotype-to-phenotype prediction using SNP data from VCF files
2. **GenoBridge_GWAS_ML** — Hybrid GWAS + ML framework for identifying candidate genes associated with complex traits

Both tools are distributed as **standalone binaries** — no Python, no dependencies, no installation required.

---

## Tools Included

| Tool | Description |
|------|-------------|
| `GenoBridge` | Ensemble ML pipeline (Ridge regression, Random Forest, XGBoost, a multi-layer perceptron, and a stacking ensemble) for phenotype prediction from genotype data |
| `GenoBridge_GWAS_ML` | Hybrid GWAS + ML pipeline for candidate gene discovery with population structure correction |

---

## Download

| Platform | GenoBridge | GenoBridge_GWAS_ML |
|----------|------------|-------------------|
| **CentOS/RHEL** | [GenoBridge_CentOS.tar.gz](https://github.com/nps-genomics/GenoBridge/releases/download/v1.0/GenoBridge_CentOS.tar.gz) | [GenoBridge_GWAS_ML_CentOS.tar.gz](https://github.com/nps-genomics/GenoBridge/releases/download/v1.0/GenoBridge_GWAS_ML_CentOS.tar.gz) |
| **Linux (Ubuntu/Debian)** | [GenoBridge_Linux.tar.gz](https://github.com/nps-genomics/GenoBridge/releases/download/v1.0/GenoBridge_Linux.tar.gz) | [GenoBridge_GWAS_ML_Linux.tar.gz](https://github.com/nps-genomics/GenoBridge/releases/download/v1.0/GenoBridge_GWAS_ML_Linux.tar.gz) |
| **macOS** | [GenoBridge_macOS.tar.gz](https://github.com/nps-genomics/GenoBridge/releases/download/v1.0/GenoBridge_macOS.tar.gz) | [GenoBridge_GWAS_ML_macOS.tar.gz](https://github.com/nps-genomics/GenoBridge/releases/download/v1.0/GenoBridge_GWAS_ML_macOS.tar.gz) |



---

## Installation

### Linux (Ubuntu/Debian)
```bash
wget https://github.com/nps-genomics/GenoBridge/releases/download/v1.0/GenoBridge_Linux.tar.gz
wget https://github.com/nps-genomics/GenoBridge/releases/download/v1.0/GenoBridge_GWAS_ML_Linux.tar.gz
tar -xzf GenoBridge_Linux.tar.gz
tar -xzf GenoBridge_GWAS_ML_Linux.tar.gz
chmod +x GenoBridge_Linux GenoBridge_GWAS_ML_Linux
sudo mv GenoBridge_Linux /usr/local/bin/GenoBridge
sudo mv GenoBridge_GWAS_ML_Linux /usr/local/bin/GenoBridge_GWAS_ML
```

### CentOS / RHEL
```bash
wget https://github.com/nps-genomics/GenoBridge/releases/download/v1.0/GenoBridge_CentOS.tar.gz
wget https://github.com/nps-genomics/GenoBridge/releases/download/v1.0/GenoBridge_GWAS_ML_CentOS.tar.gz
tar -xzf GenoBridge_CentOS.tar.gz
tar -xzf GenoBridge_GWAS_ML_CentOS.tar.gz
chmod +x GenoBridge_CentOS GenoBridge_GWAS_ML_CentOS
```

### macOS
```bash
curl -LO https://github.com/nps-genomics/GenoBridge/releases/download/v1.0/GenoBridge_macOS.tar.gz
curl -LO https://github.com/nps-genomics/GenoBridge/releases/download/v1.0/GenoBridge_GWAS_ML_macOS.tar.gz
tar -xzf GenoBridge_macOS.tar.gz
tar -xzf GenoBridge_GWAS_ML_macOS.tar.gz
chmod +x GenoBridge_macOS GenoBridge_GWAS_ML_macOS
xattr -d com.apple.quarantine GenoBridge_macOS
xattr -d com.apple.quarantine GenoBridge_GWAS_ML_macOS
```

---

## Usage

### GenoBridge
```bash
./GenoBridge_Linux \
  --phenotype data/phenotypes.csv \
  --vcf data/genotypes.vcf.gz \
  --output results/ 
  --accession-col accession_id
```

### GenoBridge_GWAS_ML
```bash
./GenoBridge_GWAS_ML_Linux \
  --phenotype data/phenotypes.csv \
  --vcf data/genotypes.vcf.gz \
  --ml-results results/phenotype_prediction_results.csv \
  --gff data/genome.gff3 \
  --output results/gwasml \
  --ml-threshold 0.3 \
  --n-pop-pcs 5
```

### Full Pipeline
```bash
# Step 1 — Genotype-to-Phenotype Prediction
./GenoBridge_Linux --phenotype phenotypes.csv --vcf genotypes.vcf.gz --output step1_output

# Step 2 — GWAS + ML Candidate Gene Discovery
./GenoBridge_GWAS_ML_Linux --phenotype phenotypes.csv --vcf genotypes.vcf.gz \
  --ml-results step1_output/phenotype_prediction_results.csv --gff genome.gff3 --output step2_output
```

---

## Input File Formats

### Phenotype CSV
accession_id,trait1,trait2,trait3
```bash
ACC001,12.5,3.2,0.87
ACC002,14.1,2.9,0.91
```


### VCF File
Standard VCF or bgzipped VCF.gz format. Sample IDs must match phenotype accession IDs.

### GFF3 (optional)
Standard GFF3 genome annotation file for candidate gene identification.

---

## Output Files

### GenoBridge

output/
├── predictions/
├── feature_importance/
├── model_performance/
└── summary_report.csv


### GenoBridge_GWAS_ML
output/
├── gwas_results/
├── candidate_genes/
├── cross_reference/
└── summary_report.csv


---

## System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| CPU | 4 cores | 16+ cores |
| RAM | 8 GB | 32–64 GB |
| Storage | 10 GB | 50 GB |
| OS | CentOS 7+, Ubuntu 18+, macOS 11+ | CentOS 9, Ubuntu 22+, macOS 13+ |

> **HPC users:** Both binaries are self-contained — no module loads or conda environments needed.

---

## Citation
Singh, N. Mendu L, Mendu V, (2026). GenoBridge: Adaptive ensemble ML pipeline for
genotype-to-phenotype prediction with integrated GWAS and candidate
gene discovery. GitHub. https://github.com/nps-genomics/GenoBridge

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

---

## Contact

**Nagendra Singh & Venugopal mendu**   
Texas A&M University 
Email: nagendra.singh@agnet.tamu.edu, Venugopal.Mendu@tamuk.edu
GitHub: [@nps-genomics](https://github.com/nps-genomics)
