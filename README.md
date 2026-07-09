# GenoBridge

**An automated, sample-size-adaptive machine-learning pipeline for genomic
prediction and mixed-model (GWAS) association, with a predictability gate that
links the two.**

GenoBridge takes standard genotype (VCF) and phenotype (CSV) inputs and, without
manual tuning: (1) selects prediction models appropriate to the dataset size,
(2) evaluates every trait with multiple complementary metrics, (3) runs a
kinship linear mixed-model GWAS only for traits that show heritable signal, and
(4) annotates the resulting associations to nearby genes.

- **Open source** (Apache 2.0) — full Python source, no hidden binaries.
- **Reproducible** — one-command conda environment; a worked example is included.
- **Statistically grounded** — reports Pearson *r*, Spearman ρ, RMSE, and
  cross-validated R²; GWAS uses a kinship mixed model with per-trait genomic
  inflation (λ) reporting.

---

## Contents
- [How it works](#how-it-works)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Tool 1 — genomic prediction](#tool-1--genomic-prediction-genobridge_fixedpy)
- [Tool 2 — mixed-model GWAS](#tool-2--mixed-model-gwas-genobridge_gwas_gemmapy)
- [Input formats](#input-formats)
- [Output files](#output-files)
- [Worked example](#worked-example)
- [Dependencies](#dependencies)
- [License and citation](#license-and-citation)
- [FAQ / troubleshooting](#faq--troubleshooting)

---

## How it works

GenoBridge runs in two steps that share one idea: **spend GWAS effort only on
traits the genome can actually predict.**

```
  phenotypes.csv ┐
                 ├─► [1] PREDICTION ─► per-trait accuracy (r, ρ, RMSE, CV-R²)
  genotypes.vcf ─┘         │
                           ▼
                     predictability gate  (best CV-selected r ≥ 0.30 ?)
                           │ yes
                           ▼
  genotypes.vcf ─► [2] MIXED-MODEL GWAS (GEMMA) ─► associations + gene annotation
```

**Step 1** classifies the dataset into a size regime and, for each trait, trains
five models (Ridge, Random Forest, XGBoost, a neural network, and a stacking
ensemble), selecting the best per trait by **cross-validated R²** (not test-set
performance, to avoid selection bias). It reports Pearson *r*, Spearman ρ, RMSE,
and CV-R² for every trait and model.

**Step 2** runs a **kinship linear mixed-model GWAS** (via GEMMA) for each trait
that passes the gate, controlling for population structure, and annotates
significant SNPs to the nearest gene. Genomic inflation (λ) is reported per
trait so structure control is transparent.

> The gate is a **heritability screen**, not a claim about which loci matter:
> a trait with negligible genome-wide predictability is unlikely to reward
> association testing, so GWAS effort is focused where signal exists. Prediction
> accuracy and GWAS detectability are distinct — a trait may have a strong
> single-locus association yet modest genome-wide predictability, or vice versa.

---

## Installation

GenoBridge is Python, and its GWAS module calls two external command-line tools
(**GEMMA** and **PLINK2**). The easiest setup installs everything at once:

```bash
git clone https://github.com/nps-genomics/GenoBridge.git
cd GenoBridge
conda env create -f environment.yml
conda activate genobridge
```

Verify the external tools are visible:

```bash
gemma  -h    # should print GEMMA version >= 0.98.5
plink2 --version
```

<details>
<summary>Manual / pip installation</summary>

```bash
pip install -r requirements.txt          # Python packages only
conda install -c bioconda gemma plink2   # external tools (not pip-installable)
```
Python packages: numpy, pandas, scipy, scikit-learn, xgboost, matplotlib.
</details>

---

## Quick start

```bash
# Step 1 — prediction + trait triage
python genobridge_FIXED.py \
  --phenotype examples/rice_phenotype.csv \
  --vcf       examples/rice_geno.vcf \
  --accession-col accession_id \
  --output    results/ \
  --no-outlier-removal

# Step 2 — mixed-model GWAS on the traits that passed the gate
python genobridge_gwas_gemma.py \
  --label     rice \
  --phenotype examples/rice_phenotype.csv \
  --vcf       examples/rice_geno.vcf \
  --gff       examples/rice_genes.gff3 \
  --ml-results results/phenotype_prediction_results.csv \
  --output    results/gwas/
```

Results appear in `results/` (prediction) and `results/gwas/` (associations,
figures, and per-trait significant-hit tables).

---

## Tool 1 — genomic prediction (`genobridge_FIXED.py`)

Trains and evaluates prediction models per trait and writes a per-trait accuracy
table used as the GWAS gate.

```
python genobridge_FIXED.py --phenotype FILE --vcf FILE [options]
```

| Argument | Default | Description |
|---|---|---|
| `--phenotype` | *(required)* | Phenotype CSV (accession column + trait columns) |
| `--vcf` | *(required)* | Genotype VCF |
| `--accession-col` | `accession_id` | Name of the accession-ID column |
| `--exclude-cols` | `replicate_id` | Non-trait columns to ignore |
| `--output` | `results` | Output directory |
| `--no-outlier-removal` | off | **Recommended for multi-trait data.** Disables per-trait |z|>3 outlier removal, which can drop many accessions across many traits |
| `--max-missing-frac` | `0.30` | Traits with more missingness than this are dropped and reported, not imputed (imputing high-missing traits fabricates data) |
| `--n-components` | auto | PCA components (default set by regime) |
| `--test-size` | auto | Held-out test fraction (default set by regime) |
| `--force-regime` | *(unset)* | **Advanced/analysis only.** Overrides automatic regime selection; normal users leave this unset |

> **Tip:** for datasets with many traits (e.g. multi-environment yield trials),
> use `--no-outlier-removal`. Per-trait outlier filtering across dozens of
> heterogeneous traits can remove a large fraction of accessions.

---

## Tool 2 — mixed-model GWAS (`genobridge_gwas_gemma.py`)

Runs a kinship linear mixed-model GWAS (GEMMA) for gated traits and annotates
hits to nearby genes.

```
python genobridge_gwas_gemma.py --phenotype FILE --vcf FILE --gff FILE [options]
```

| Argument | Default | Description |
|---|---|---|
| `--phenotype` | *(required)* | Phenotype CSV (same format as Tool 1) |
| `--vcf` | *(required)* | Genotype VCF |
| `--gff` | *(required)* | Gene annotation (GFF3) for nearest-gene labeling |
| `--ml-results` | *(none)* | Prediction results CSV from Tool 1; enables the gate |
| `--gate` | `0.30` | Predictability threshold; traits with best *r* ≥ this go to GWAS. Set `0` to run all traits |
| `--gene-function` | *(none)* | Optional `gene_id,function` CSV to add/override gene functions (GFF `Note=`/`description=` used automatically otherwise) |
| `--accession-col` | `accession_id` | Accession-ID column name |
| `--maf` | `0.05` | Minor-allele-frequency filter |
| `--geno-missing` | `0.10` | Per-SNP missingness filter |
| `--label` | `gwas` | Cosmetic label for the report/figures |
| `--output` | `gwas_gemma_out` | Output directory |

**What it does per run:** converts the VCF to PLINK binary (MAF/missingness
filtered), builds a centered genomic relatedness matrix, fits a mixed-model
association test per gated trait, annotates significant SNPs to the nearest
gene, and reports genomic inflation (λ) for each trait.

---

## Input formats

**Phenotype CSV** — one accession-ID column plus one column per trait. Missing
values blank or `NA`.

```
accession_id,Brix,Amylose,Plant_height
S001,18.2,22.1,95
S002,17.5,,102
S003,,20.8,88
```

**Genotype VCF** — standard VCF (bgzipped or plain). Sample IDs in the VCF header
must match the `accession_id` values in the phenotype CSV. Multi-allelic sites
are filtered automatically; only biallelic SNPs are used.

**GFF3** (GWAS only) — standard gene annotation. Functional descriptions are read
from `Note=`/`description=`/`product=` attributes when present.

**Optional gene-function CSV** — `gene_id,function`, to supply or override
functional labels for annotated genes.

---

## Output files

**Prediction (`results/`)**
- `phenotype_prediction_results.csv` — per trait: Pearson *r*, Spearman ρ, RMSE,
  CV-R² for each model, plus the CV-selected best model and its metrics. This is
  the file passed to the GWAS gate.
- `predictions.csv` — observed vs predicted values on the held-out test set.

**GWAS (`results/gwas/`)**
- `gwas_summary.csv` — per gated trait: λ, and Bonferroni/suggestive hit counts.
- `assoc/<trait>.sig.csv` — significant SNPs with clear column names, sorted by
  significance: `Chromosome, SNP_ID, Position_bp, Effect_allele, Ref_allele,
  Allele_freq, Effect_size, SE, P_value, Neg_log10_P, Lambda_REML, Nearest_gene,
  Distance_bp, Gene_function`.
- `assoc/<trait>.assoc.txt` — full raw GEMMA output (complete record).
- `figures/` — Manhattan and QQ plots per trait (QQ titles show λ), plus a
  per-trait λ comparison figure.

---

## Worked example

The `examples/` directory contains a small demonstration dataset. From the repo
root:

```bash
# prediction
python genobridge_FIXED.py \
  --phenotype examples/rice_phenotype.csv \
  --vcf examples/rice_geno.vcf \
  --output example_results/ --no-outlier-removal

# GWAS on gated traits
python genobridge_gwas_gemma.py --label example \
  --phenotype examples/rice_phenotype.csv \
  --vcf examples/rice_geno.vcf \
  --gff examples/rice_genes.gff3 \
  --ml-results example_results/phenotype_prediction_results.csv \
  --output example_results/gwas/
```

Expected: a prediction table with all four metrics, and — for traits that pass
the gate — Manhattan/QQ plots with λ near 1.0 and annotated significant SNPs.

---

## Dependencies

- **Python ≥ 3.8** with numpy, pandas, scipy, scikit-learn, xgboost, matplotlib
- **GEMMA ≥ 0.98.5** and **PLINK2** on the system PATH (GWAS module)

All installable with `conda env create -f environment.yml`.

---

## License and citation

GenoBridge is released under the **Apache License 2.0** (see `LICENSE`). You may
use, modify, and redistribute it provided you retain the attribution notices
(`NOTICE`) and, for academic use, cite the publication below. Modified versions
must state that they were changed and preserve the original attribution.

**If you use GenoBridge in academic work, please cite:**

> Singh, N.P., Mendu, V. (2026). *GenoBridge: an automated
> sample-size-adaptive machine-learning pipeline for genomic prediction and
> ML-gated GWAS.* Plant Methods. DOI: *[waiting on acceptance]*

and the archived release (Zenodo DOI: *[waiting on acceptance]*). A
machine-readable citation is in `CITATION.cff`.

---

## FAQ / troubleshooting

**"gemma: command not found" / "plink2: command not found"** — the external
tools aren't on your PATH. Install with `conda install -c bioconda gemma plink2`
and re-activate the environment.

**Wheat-scale / large VCFs run out of memory in the prediction step** — the
prediction pipeline loads genotypes into memory. For very large panels,
LD-prune the VCF first (e.g. `plink2 --indep-pairwise 100 10 0.2`) and run
prediction on the pruned set; the GWAS step streams the full VCF and is
unaffected.

**All my GWAS p-values look null / λ ≈ 0** — ensure the output directory is
fresh (delete stale `geno.*` files) and that phenotype accession IDs exactly
match the VCF sample IDs. A near-zero overlap silently produces null results.

**Which traits went to GWAS?** — only those with best CV-selected *r* ≥ `--gate`
(default 0.30). Use `--gate 0` to run all traits regardless.

---

*.*
