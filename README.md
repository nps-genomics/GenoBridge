# GenoBridge
Adaptive ensemble ML pipeline for genotype-to-phenotype prediction  with integrated GWAS and candidate gene discovery. Supports 30->100,000+  samples 

# AdaptGenoPred
<<<<<<< Updated upstream
Adaptive ensemble ML pipeline for genotype-to-phenotype prediction  with integrated GWAS and candidate gene discovery. Supports 30->100,000+  samples 
=======

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20CentOS%2FRHEL-blue)
![Version](https://img.shields.io/badge/version-1.0-green)

**Adaptive ensemble ML pipeline for genotype-to-phenotype prediction with integrated GWAS and candidate gene discovery.**  
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
- [Examples](#examples)
- [System Requirements](#system-requirements)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Overview

AdaptGenoPred is a two-tool pipeline designed for plant genomics research:

1. **AdaptGenoPred** — Adaptive ensemble machine learning for genotype-to-phenotype prediction using SNP data from VCF files
2. **GWASml** — Hybrid GWAS + ML framework for identifying candidate genes associated with complex traits

Both tools are distributed as **standalone binaries** — no Python, no dependencies, no installation required.

---

## Tools Included

| Tool | Description |
|------|-------------|
| \`AdaptGenoPred\` | Ensemble ML pipeline (Random Forest, XGBoost, SVM) for phenotype prediction from genotype data |
| \`GWASml\` | Hybrid GWAS + ML pipeline for candidate gene discovery with population structure correction |

---

## Download

| Platform | AdaptGenoPred | GWASml |
|----------|--------------|--------|
| **CentOS/RHEL** | [AdaptGenoPred-CentOS.tar.gz](https://github.com/nps-genomics/AdaptGenoPred/releases/latest/download/AdaptGenoPred-CentOS.tar.gz) | [GWASml-CentOS.tar.gz](https://github.com/nps-genomics/AdaptGenoPred/releases/latest/download/GWASml-CentOS.tar.gz) |
| **Linux (Ubuntu/Debian)** | [AdaptGenoPred-Linux.tar.gz](https://github.com/nps-genomics/AdaptGenoPred/releases/latest/download/AdaptGenoPred-Linux.tar.gz) | [GWASml-Linux.tar.gz](https://github.com/nps-genomics/AdaptGenoPred/releases/latest/download/GWASml-Linux.tar.gz) |
| **macOS** | [AdaptGenoPred-macOS.tar.gz](https://github.com/nps-genomics/AdaptGenoPred/releases/latest/download/AdaptGenoPred-macOS.tar.gz) | [GWASml-macOS.tar.gz](https://github.com/nps-genomics/AdaptGenoPred/releases/latest/download/GWASml-macOS.tar.gz) |

### SHA256 Checksums (v1.0)
\`\`\`
AdaptGenoPred-CentOS.tar.gz  1a91ae720f53bf7f3251f686b1990d4cfd51935e90b5cfb9b8faea3305fb9463
AdaptGenoPred-Linux.tar.gz   ccfe1cff93a71598fc453483e0bc5c37bba7e04fe72c2c86605d31daeac4af37
AdaptGenoPred-macOS.tar.gz   255557600a5595e36900c6142b426d0127980099fe29acc1246784505cf8a741
GWASml-CentOS.tar.gz         5f495ccb632a6f126c055bb3ed93245f33427d66516f872ba6f0dd58cdc8219d
GWASml-Linux.tar.gz          a5a195ca350f9bd184b8867907f1bd43a5c9540871b49fea14cc6ff80acaee4c
GWASml-macOS.tar.gz          a96d61efd539a8837f19a767e0361d981a5a5359501b2df74c9fda7109856480
\`\`\`

---

## Installation

### Linux (Ubuntu/Debian)
\`\`\`bash
wget https://github.com/nps-genomics/AdaptGenoPred/releases/latest/download/AdaptGenoPred-Linux.tar.gz
wget https://github.com/nps-genomics/AdaptGenoPred/releases/latest/download/GWASml-Linux.tar.gz
tar -xzf AdaptGenoPred-Linux.tar.gz
tar -xzf GWASml-Linux.tar.gz
chmod +x AdaptGenoPred-Linux GWASml-Linux
./AdaptGenoPred-Linux 
./GWASml-Linux
or
sudo mv AdaptGenoPred-Linux /usr/local/bin/AdaptGenoPred
sudo mv GWASml-Linux /usr/local/bin/GWASml
\`\`\`

### CentOS / RHEL
\`\`\`bash
wget https://github.com/nps-genomics/AdaptGenoPred/releases/latest/download/AdaptGenoPred-CentOS.tar.gz
wget https://github.com/nps-genomics/AdaptGenoPred/releases/latest/download/GWASml-CentOS.tar.gz
tar -xzf AdaptGenoPred-CentOS.tar.gz
tar -xzf GWASml-CentOS.tar.gz
chmod +x AdaptGenoPred-CentOS GWASml-CentOS
\`\`\`

### macOS
\`\`\`bash
curl -LO https://github.com/nps-genomics/AdaptGenoPred/releases/latest/download/AdaptGenoPred-macOS.tar.gz
curl -LO https://github.com/nps-genomics/AdaptGenoPred/releases/latest/download/GWASml-macOS.tar.gz
tar -xzf AdaptGenoPred-macOS.tar.gz
tar -xzf GWASml-macOS.tar.gz
chmod +x AdaptGenoPred-macOS GWASml-macOS
xattr -d com.apple.quarantine AdaptGenoPred-macOS
xattr -d com.apple.quarantine GWASml-macOS
\`\`\`

---

## Usage

### AdaptGenoPred
\`\`\`bash
./AdaptGenoPred-Linux \
  --phenotype data/phenotypes.csv \
  --vcf data/genotypes.vcf.gz \
  --output results/ \
  --n-pop-pcs 5 \
  --accession-col accession_id
\`\`\`

### GWASml
\`\`\`bash
./GWASml-Linux \
  --phenotype data/phenotypes.csv \
  --vcf data/genotypes.vcf.gz \
  --ml-results results/summary_report.csv \
  --gff data/genome.gff3 \
  --output results/gwasml \
  --ml-threshold 0.3 \
  --n-pop-pcs 5
\`\`\`

### Full Pipeline
\`\`\`bash
# Step 1
./AdaptGenoPred-Linux --phenotype phenotypes.csv --vcf genotypes.vcf.gz --output step1_output

# Step 2
./GWASml-Linux --phenotype phenotypes.csv --vcf genotypes.vcf.gz \
  --ml-results step1_output/summary_report.csv --gff genome.gff3 --output step2_output
\`\`\`

---

## Input File Formats

### Phenotype CSV
\`\`\`
accession_id,trait1,trait2,trait3
ACC001,12.5,3.2,0.87
ACC002,14.1,2.9,0.91
\`\`\`

### VCF File
Standard VCF or bgzipped VCF.gz format. Sample IDs must match phenotype accession IDs.

### GFF3 (optional)
Standard GFF3 genome annotation file for candidate gene identification.

---

## Output Files

### AdaptGenoPred
\`\`\`
output/
├── predictions/
├── feature_importance/
├── model_performance/
└── summary_report.csv
\`\`\`

### GWASml
\`\`\`
output/
├── gwas_results/
├── candidate_genes/
├── cross_reference/
└── summary_report.csv
\`\`\`

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

\`\`\`
Singh, N. (2026). AdaptGenoPred: Adaptive ensemble ML pipeline for
genotype-to-phenotype prediction with integrated GWAS and candidate
gene discovery. GitHub. https://github.com/nps-genomics/AdaptGenoPred
\`\`\`

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

---

## Contact

**Nagendra Singh**  
Caesar Kleberg Wildlife Research Institute  
Texas A&M University Kingsville  
GitHub: [@nps-genomics](https://github.com/nps-genomics)
>>>>>>> Stashed changes
