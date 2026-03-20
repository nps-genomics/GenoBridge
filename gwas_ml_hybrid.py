"""
GWAS-ML Hybrid Pipeline for SNP-Trait Association Analysis.

Performs genome-wide association analysis, cross-references significant loci
with ML prediction performance, and annotates candidate genes using GFF3
annotations (TAIR10, MSU7, or any species).

Usage (command line):
    python gwas_ml_hybrid.py \\
        --phenotype data/phenotypes.csv \\
        --vcf data/genotypes.vcf \\
        --ml-results results/phenotype_prediction_results.csv \\
        --output gwas_results/

    python gwas_ml_hybrid.py \\
        --phenotype data/rice_traits.csv \\
        --vcf data/rice_44k.vcf \\
        --ml-results results_rice/phenotype_prediction_results.csv \\
        --gff data/rice_MSU7_genes.gff3 \\
        --accession-col HybID \\
        --n-pop-pcs 5 \\
        --ml-threshold 0.3 \\
        --output gwas_rice/

    For full options:
        python gwas_ml_hybrid.py --help

Usage (Python API):
    from gwas_ml_hybrid import GWASMLHybrid

    hybrid = GWASMLHybrid(
        phenotype_file="pheno.csv",
        vcf_file="geno.vcf",
        ml_results_file="results/phenotype_prediction_results.csv",
        gff_file="genes.gff3",
    )
    hybrid.run()

Requirements:
    pip install numpy pandas scipy scikit-learn scikit-allel matplotlib

Output (in output directory):
    - Per-trait GWAS results CSVs with gene annotations
    - Manhattan plots and QQ plots
    - Summary table crossing ML performance with GWAS hits
    - Highlight table for manuscript (top traits + top SNPs + genes)
"""

import pandas as pd
import numpy as np
import allel
import os
import sys
import gzip
import shutil
import logging
import warnings
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler
from sklearn.impute import KNNImputer
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict
from typing import Optional

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


# ═══════════════════════════════════════════════════════════════════════════
# 1. Data loading (mirrors the adaptive pipeline logic)
# ═══════════════════════════════════════════════════════════════════════════

def load_phenotype(phenotype_file: str, accession_col: str = "accession_id",
                   exclude_cols: list = None, remove_outliers: bool = True) -> tuple:
    """Load and process phenotype data (same as adaptive pipeline)."""
    exclude_cols = exclude_cols or ["replicate_id"]
    df = pd.read_csv(phenotype_file)
    skip = set([accession_col] + exclude_cols)
    traits = [c for c in df.columns if c not in skip]

    df_agg = df.groupby(accession_col)[traits].mean().reset_index()

    # KNN imputation
    imp = KNNImputer(n_neighbors=min(5, max(2, len(df_agg) // 20)))
    df_agg[traits] = imp.fit_transform(df_agg[traits])

    # Robust scaling
    scaler = RobustScaler()
    df_agg[traits] = scaler.fit_transform(df_agg[traits])

    # Outlier removal
    if remove_outliers:
        before = len(df_agg)
        for trait in traits:
            mu, sigma = df_agg[trait].mean(), df_agg[trait].std()
            if sigma > 0:
                df_agg = df_agg[np.abs((df_agg[trait] - mu) / sigma) <= 3]
        logging.info(f"Phenotype: {len(df_agg)} accessions ({before}→{len(df_agg)} after outlier removal), {len(traits)} traits")
    else:
        logging.info(f"Phenotype: {len(df_agg)} accessions (outlier removal skipped), {len(traits)} traits")

    return df_agg, traits


def load_genotype(vcf_file: str, pheno_ids: list) -> tuple:
    """Load genotype data and match accessions."""
    logging.info(f"Loading VCF: {vcf_file}")
    callset = allel.read_vcf(
        vcf_file,
        fields=["samples", "calldata/GT", "variants/ID", "variants/CHROM", "variants/POS"],
    )
    if callset is None:
        raise ValueError("Failed to read VCF")

    genotypes = callset["calldata/GT"].sum(axis=2).astype(float).T
    genotypes[genotypes < 0] = np.nan
    sample_ids = [str(s) for s in callset["samples"]]
    snp_ids = list(callset["variants/ID"])
    chroms = list(callset["variants/CHROM"])
    positions = list(callset["variants/POS"])

    logging.info(f"VCF: {len(sample_ids)} samples, {len(snp_ids)} SNPs")

    # Handle duplicate SNP IDs
    if len(snp_ids) != len(set(snp_ids)):
        logging.warning("Duplicate SNP IDs — generating unique IDs")
        snp_ids = [f"SNP_{i}" for i in range(len(snp_ids))]

    # Build SNP info dataframe
    snp_info = pd.DataFrame({
        "SNP_ID": snp_ids,
        "CHROM": chroms,
        "POS": positions,
    })

    # Match accessions (same multi-strategy as adaptive pipeline)
    pheno_set = set(str(a) for a in pheno_ids)
    vcf_set = set(sample_ids)
    direct = pheno_set & vcf_set
    logging.info(f"ID matching — direct: {len(direct)} matches")

    if len(direct) >= max(10, len(pheno_set) * 0.3):
        mask = [s in pheno_set for s in sample_ids]
        matched_ids = [s for s in sample_ids if s in pheno_set]
    else:
        # Numeric matching
        vcf_num = {}
        for i, vid in enumerate(sample_ids):
            try:
                vcf_num[str(int(vid))] = i
            except ValueError:
                continue
        pheno_num = {str(int(p)): p for p in pheno_set if p.isdigit()}
        num_matches = set(vcf_num) & set(pheno_num)
        logging.info(f"ID matching — numeric: {len(num_matches)} matches")

        if len(num_matches) >= 10:
            indices = [vcf_num[k] for k in num_matches]
            matched_ids = [pheno_num[k] for k in num_matches]
            mask = [False] * len(sample_ids)
            for idx in indices:
                mask[idx] = True
        else:
            raise ValueError("Could not match accession IDs between VCF and phenotype data.")

    geno_matrix = genotypes[mask]

    # Impute missing genotype calls
    imp = KNNImputer(n_neighbors=min(5, max(2, len(matched_ids) // 20)))
    geno_matrix = imp.fit_transform(geno_matrix)

    logging.info(f"Matched: {geno_matrix.shape[0]} accessions, {geno_matrix.shape[1]} SNPs")
    return geno_matrix, matched_ids, snp_info


# ═══════════════════════════════════════════════════════════════════════════
# 2. GWAS — per-SNP linear regression with population structure correction
# ═══════════════════════════════════════════════════════════════════════════

def compute_population_pcs(geno_matrix: np.ndarray, n_pcs: int = 5) -> np.ndarray:
    """Compute top PCs from genotype matrix for population structure correction."""
    logging.info(f"Computing {n_pcs} population structure PCs")
    pca = PCA(n_components=n_pcs)
    pcs = pca.fit_transform(geno_matrix)
    var_explained = pca.explained_variance_ratio_.sum()
    logging.info(f"Top {n_pcs} PCs explain {var_explained:.1%} of genotype variance")
    return pcs


def run_gwas_single_trait(y: np.ndarray, geno_matrix: np.ndarray,
                          covariates: np.ndarray) -> tuple:
    """Run GWAS for one trait: linear regression per SNP with covariates.

    Model: y = beta_0 + beta_snp * SNP_j + beta_cov * covariates + epsilon

    Uses vectorized matrix operations for speed.

    Returns: (betas, standard_errors, p_values) arrays of shape (n_snps,)
    """
    n = len(y)
    n_snps = geno_matrix.shape[1]

    # Prepare covariate matrix (intercept + population PCs)
    C = np.column_stack([np.ones(n), covariates])  # (n, 1+n_pcs)
    n_cov = C.shape[1]

    # Project out covariates from y
    # y_adj = y - C @ (C^T C)^{-1} C^T y
    CtC_inv = np.linalg.inv(C.T @ C)
    hat_matrix = C @ CtC_inv @ C.T
    y_adj = y - hat_matrix @ y
    y_ss = y_adj.T @ y_adj  # sum of squares of adjusted y

    # For each SNP, project out covariates from SNP, then regress
    # This is equivalent to the partial regression coefficient
    betas = np.zeros(n_snps)
    se = np.zeros(n_snps)
    pvals = np.ones(n_snps)

    # Process in chunks for memory efficiency
    chunk_size = 10000
    df_resid = n - n_cov - 1  # degrees of freedom

    if df_resid <= 0:
        logging.warning(f"Not enough degrees of freedom (n={n}, cov={n_cov})")
        return betas, se, pvals

    for start in range(0, n_snps, chunk_size):
        end = min(start + chunk_size, n_snps)
        G_chunk = geno_matrix[:, start:end]  # (n, chunk)

        # Adjust SNPs for covariates
        G_adj = G_chunk - hat_matrix @ G_chunk  # (n, chunk)

        # Variance of adjusted SNPs
        g_ss = np.sum(G_adj ** 2, axis=0)  # (chunk,)

        # Avoid division by zero for invariant SNPs
        valid = g_ss > 1e-10

        # Beta = (G_adj^T y_adj) / (G_adj^T G_adj)
        numerator = G_adj.T @ y_adj  # (chunk,)
        b = np.zeros(end - start)
        b[valid] = numerator[valid] / g_ss[valid]
        betas[start:end] = b

        # Residual sum of squares
        rss = np.zeros(end - start)
        rss[valid] = y_ss - (numerator[valid] ** 2) / g_ss[valid]
        rss[~valid] = y_ss

        # Standard error
        mse = rss / df_resid
        s = np.zeros(end - start)
        s[valid] = np.sqrt(np.maximum(mse[valid] / g_ss[valid], 1e-30))
        se[start:end] = s

        # t-statistic and p-value
        t_stat = np.zeros(end - start)
        t_stat[valid] = b[valid] / np.maximum(s[valid], 1e-30)
        p = np.ones(end - start)
        p[valid] = 2 * stats.t.sf(np.abs(t_stat[valid]), df_resid)
        pvals[start:end] = p

    return betas, se, pvals


def run_gwas_all_traits(pheno_df: pd.DataFrame, traits: list,
                        geno_matrix: np.ndarray, matched_ids: list,
                        snp_info: pd.DataFrame, pop_pcs: np.ndarray,
                        output_dir: str) -> pd.DataFrame:
    """Run GWAS for all traits and save results."""
    os.makedirs(output_dir, exist_ok=True)

    # Align phenotype data with matched genotype order
    pheno_df = pheno_df.copy()
    pheno_df["accession_id"] = pheno_df["accession_id"].astype(str)
    pheno_aligned = pheno_df.set_index("accession_id").loc[matched_ids]

    n_snps = geno_matrix.shape[1]
    bonferroni = 0.05 / n_snps
    suggestive = 1e-5

    logging.info(f"Running GWAS: {len(traits)} traits × {n_snps} SNPs")
    logging.info(f"Bonferroni threshold: p < {bonferroni:.2e}")
    logging.info(f"Suggestive threshold: p < {suggestive:.2e}")

    summary_rows = []

    for i, trait in enumerate(traits):
        y = pheno_aligned[trait].values.astype(float)

        # Skip if no variance
        if np.std(y) < 1e-10:
            logging.warning(f"Skipping {trait}: no variance")
            continue

        betas, se, pvals = run_gwas_single_trait(y, geno_matrix, pop_pcs)

        # Count hits
        n_bonf = np.sum(pvals < bonferroni)
        n_sugg = np.sum(pvals < suggestive)

        # Save per-trait results (only suggestive and above to save space)
        mask_save = pvals < suggestive
        if mask_save.sum() > 0:
            trait_df = snp_info.loc[mask_save].copy()
            trait_df["Beta"] = betas[mask_save]
            trait_df["SE"] = se[mask_save]
            trait_df["P_value"] = pvals[mask_save]
            trait_df["-log10P"] = -np.log10(np.maximum(pvals[mask_save], 1e-300))
            trait_df["Trait"] = trait
            trait_df = trait_df.sort_values("P_value")
            trait_df.to_csv(
                os.path.join(output_dir, f"gwas_{trait}.csv"), index=False
            )
        else:
            trait_df = pd.DataFrame()

        # Save ALL p-values for Manhattan plot (compressed)
        all_pvals_df = pd.DataFrame({
            "CHROM": snp_info["CHROM"].values,
            "POS": snp_info["POS"].values,
            "P": pvals,
        })
        all_pvals_df.to_csv(
            os.path.join(output_dir, f"gwas_pvals_{trait}.csv.gz"),
            index=False, compression="gzip",
        )

        min_p = pvals.min()
        summary_rows.append({
            "Trait": trait,
            "N_Bonferroni": n_bonf,
            "N_Suggestive": n_sugg,
            "Min_P": min_p,
            "Min_P_log10": -np.log10(max(min_p, 1e-300)),
            "Top_SNP_idx": np.argmin(pvals),
        })

        if (i + 1) % 10 == 0 or i == 0:
            logging.info(f"  [{i+1}/{len(traits)}] {trait}: "
                         f"Bonferroni={n_bonf}, suggestive={n_sugg}, "
                         f"min_p={min_p:.2e}")

    summary = pd.DataFrame(summary_rows)

    # Add SNP info for top hit per trait
    summary["Top_SNP"] = summary["Top_SNP_idx"].apply(
        lambda idx: snp_info.iloc[idx]["SNP_ID"])
    summary["Top_CHROM"] = summary["Top_SNP_idx"].apply(
        lambda idx: snp_info.iloc[idx]["CHROM"])
    summary["Top_POS"] = summary["Top_SNP_idx"].apply(
        lambda idx: snp_info.iloc[idx]["POS"])
    summary = summary.drop(columns=["Top_SNP_idx"])

    summary.to_csv(os.path.join(output_dir, "gwas_summary.csv"), index=False)
    logging.info(f"GWAS complete. Summary saved.")
    return summary


# ═══════════════════════════════════════════════════════════════════════════
# 3. Gene annotation using TAIR10 GFF3
# ═══════════════════════════════════════════════════════════════════════════

def download_tair10_gff(output_path: str = "TAIR10_GFF3_genes.gff"):
    """Download TAIR10 gene annotation if not present."""
    if os.path.exists(output_path):
        logging.info(f"TAIR10 GFF found: {output_path}")
        return output_path

    logging.info("Downloading TAIR10 GFF3 annotation...")
    import urllib.request
    url = "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff"
    try:
        urllib.request.urlretrieve(url, output_path)
        logging.info(f"Downloaded to {output_path}")
    except Exception as e:
        logging.warning(f"Auto-download failed: {e}")
        logging.warning("Please download TAIR10_GFF3_genes.gff manually from "
                        "https://www.arabidopsis.org/download/list?dir=Genes/TAIR10_genome_release/TAIR10_gff3")
        logging.warning("Or provide a GFF3 file via --gff argument")
        return None
    return output_path


def load_gene_annotations(gff_file: str) -> pd.DataFrame:
    """Parse GFF3 to extract gene positions and descriptions.

    Handles TAIR10, MSU7, Grapedia (PN40024 v3), and other standard GFF3 formats.
    Captures ID, Name, Alias, and Note attributes.
    """
    if gff_file is None or not os.path.exists(gff_file):
        logging.warning("No GFF file available. Gene annotation will be skipped.")
        return None

    logging.info(f"Parsing gene annotations from {gff_file}")
    genes = []

    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "gene":
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = fields[8]

            # Parse attributes
            gene_id = ""
            gene_name = ""
            alias = ""
            description = ""
            for attr in attrs.split(";"):
                if attr.startswith("ID="):
                    gene_id = attr.split("=", 1)[1]
                elif attr.startswith("Name="):
                    gene_name = attr.split("=", 1)[1]
                elif attr.startswith("Alias="):
                    alias = attr.split("=", 1)[1]
                elif attr.startswith("Note="):
                    description = attr.split("=", 1)[1].replace("%20", " ").replace("%2C", ",")

            # Use alias as gene_name fallback if Name is absent
            if not gene_name and alias:
                gene_name = alias

            genes.append({
                "Gene_ID": gene_id,
                "Gene_Name": gene_name,
                "CHROM": chrom,
                "Gene_Start": start,
                "Gene_End": end,
                "Strand": strand,
                "Description": description,
            })

    gene_df = pd.DataFrame(genes)
    logging.info(f"Loaded {len(gene_df)} gene annotations")
    return gene_df


def annotate_snps_with_genes(snp_df: pd.DataFrame, gene_df: pd.DataFrame,
                              window: int = 10000) -> pd.DataFrame:
    """Find nearest gene(s) for each SNP within a window (default ±10kb)."""
    if gene_df is None or snp_df.empty:
        return snp_df

    # Normalize chromosome names
    snp_chrom = snp_df["CHROM"].astype(str)
    gene_chrom = gene_df["CHROM"].astype(str)

    # Robust harmonization: strip any chr/Chr/CHR prefix and leading zeros
    # to produce a bare canonical form (e.g., "chr01" -> "1", "Chr1" -> "1", "1" -> "1")
    import re

    def _normalize_chrom(s):
        """Strip chr/Chr/CHR prefix and leading zeros to get bare chromosome ID."""
        s = re.sub(r'^[Cc][Hh][Rr]', '', str(s))  # strip chr/Chr/CHR
        s = s.lstrip('0') or '0'                    # strip leading zeros (keep "0" for chr00)
        return s

    snp_df = snp_df.copy()
    gene_df = gene_df.copy()

    snp_df["CHROM_match"] = snp_df["CHROM"].astype(str).map(_normalize_chrom)
    gene_df["CHROM_match"] = gene_df["CHROM"].astype(str).map(_normalize_chrom)

    # Log the harmonization
    snp_sample = snp_chrom.iloc[0] if len(snp_chrom) > 0 else "?"
    gene_sample = gene_chrom.iloc[0] if len(gene_chrom) > 0 else "?"
    snp_norm = snp_df["CHROM_match"].iloc[0] if len(snp_df) > 0 else "?"
    gene_norm = gene_df["CHROM_match"].iloc[0] if len(gene_df) > 0 else "?"
    logging.info(
        f"Chromosome harmonization: SNP '{snp_sample}' → '{snp_norm}', "
        f"GFF '{gene_sample}' → '{gene_norm}'"
    )

    # Build gene index per chromosome
    gene_by_chrom = {}
    for chrom in gene_df["CHROM_match"].unique():
        sub = gene_df[gene_df["CHROM_match"] == chrom].copy()
        sub = sub.sort_values("Gene_Start")
        gene_by_chrom[chrom] = sub

    nearest_genes = []
    nearest_names = []
    nearest_descs = []
    distances = []
    locations = []  # upstream, downstream, within

    for _, snp in snp_df.iterrows():
        chrom = snp["CHROM_match"]
        pos = int(snp["POS"])

        if chrom not in gene_by_chrom:
            nearest_genes.append("")
            nearest_names.append("")
            nearest_descs.append("")
            distances.append(np.nan)
            locations.append("")
            continue

        genes = gene_by_chrom[chrom]

        # Find genes within window
        nearby = genes[
            (genes["Gene_Start"] - window <= pos) &
            (genes["Gene_End"] + window >= pos)
        ]

        if len(nearby) == 0:
            # Find absolute nearest
            mid = (genes["Gene_Start"] + genes["Gene_End"]) / 2
            dists = np.abs(mid - pos)
            idx = dists.idxmin()
            nearest_genes.append(genes.loc[idx, "Gene_ID"])
            nearest_names.append(genes.loc[idx, "Gene_Name"])
            nearest_descs.append(genes.loc[idx, "Description"])
            distances.append(int(dists.loc[idx]))
            locations.append("distal")
        else:
            # Pick the gene where SNP is inside, or nearest edge
            best_idx = None
            best_dist = float("inf")
            best_loc = ""

            for idx, gene in nearby.iterrows():
                if gene["Gene_Start"] <= pos <= gene["Gene_End"]:
                    best_idx = idx
                    best_dist = 0
                    best_loc = "within"
                    break
                else:
                    d = min(abs(pos - gene["Gene_Start"]),
                            abs(pos - gene["Gene_End"]))
                    if d < best_dist:
                        best_dist = d
                        best_idx = idx
                        if pos < gene["Gene_Start"]:
                            best_loc = "upstream"
                        else:
                            best_loc = "downstream"

            nearest_genes.append(genes.loc[best_idx, "Gene_ID"])
            nearest_names.append(genes.loc[best_idx, "Gene_Name"])
            nearest_descs.append(genes.loc[best_idx, "Description"])
            distances.append(int(best_dist))
            locations.append(best_loc)

    snp_df["Nearest_Gene_ID"] = nearest_genes
    snp_df["Nearest_Gene_Name"] = nearest_names
    snp_df["Gene_Description"] = nearest_descs
    snp_df["Distance_to_Gene"] = distances
    snp_df["SNP_Location"] = locations

    if "CHROM_match" in snp_df.columns:
        snp_df = snp_df.drop(columns=["CHROM_match"])

    return snp_df


# ═══════════════════════════════════════════════════════════════════════════
# 4. Cross-reference GWAS with ML results
# ═══════════════════════════════════════════════════════════════════════════

def cross_reference_gwas_ml(gwas_summary: pd.DataFrame,
                             ml_results_file: str,
                             gwas_dir: str,
                             gene_df: pd.DataFrame,
                             output_dir: str,
                             ml_threshold: float = 0.3) -> pd.DataFrame:
    """Cross-reference GWAS hits with ML prediction performance."""
    os.makedirs(output_dir, exist_ok=True)

    ml_df = pd.read_csv(ml_results_file)

    # Merge GWAS summary with ML results
    merged = gwas_summary.merge(ml_df, on="Trait", how="left")
    merged = merged.sort_values("Best_r", ascending=False)
    merged.to_csv(os.path.join(output_dir, "gwas_ml_merged.csv"), index=False)

    # Focus on traits with ML r > threshold
    high_ml = merged[merged["Best_r"] >= ml_threshold]
    logging.info(f"Traits with ML r >= {ml_threshold}: {len(high_ml)}")

    # For each high-ML trait, load and annotate GWAS hits
    highlight_rows = []

    for _, row in high_ml.iterrows():
        trait = row["Trait"]
        gwas_file = os.path.join(gwas_dir, f"gwas_{trait}.csv")

        if not os.path.exists(gwas_file):
            continue

        trait_gwas = pd.read_csv(gwas_file)
        if len(trait_gwas) == 0:
            continue

        # Annotate with genes
        trait_gwas = annotate_snps_with_genes(trait_gwas, gene_df)

        # Save annotated version
        trait_gwas.to_csv(
            os.path.join(output_dir, f"annotated_{trait}.csv"), index=False
        )

        # Add top hits to highlight table
        top_n = min(5, len(trait_gwas))
        for j in range(top_n):
            snp_row = trait_gwas.iloc[j]
            highlight_rows.append({
                "Trait": trait,
                "ML_Best_r": row["Best_r"],
                "ML_Best_Model": row["Best_Model"],
                "SNP_ID": snp_row.get("SNP_ID", ""),
                "CHROM": snp_row.get("CHROM", ""),
                "POS": snp_row.get("POS", ""),
                "P_value": snp_row.get("P_value", ""),
                "-log10P": snp_row.get("-log10P", ""),
                "Beta": snp_row.get("Beta", ""),
                "Nearest_Gene_ID": snp_row.get("Nearest_Gene_ID", ""),
                "Nearest_Gene_Name": snp_row.get("Nearest_Gene_Name", ""),
                "Gene_Description": snp_row.get("Gene_Description", ""),
                "Distance_to_Gene": snp_row.get("Distance_to_Gene", ""),
                "SNP_Location": snp_row.get("SNP_Location", ""),
            })

    highlight = pd.DataFrame(highlight_rows)
    if not highlight.empty:
        highlight = highlight.sort_values(["ML_Best_r", "P_value"],
                                           ascending=[False, True])
        highlight.to_csv(os.path.join(output_dir, "highlight_table.csv"),
                         index=False)
        logging.info(f"Highlight table: {len(highlight)} entries for "
                     f"{highlight['Trait'].nunique()} traits")
    else:
        logging.warning("No GWAS hits found for high-ML traits at suggestive threshold")

    return highlight


# ═══════════════════════════════════════════════════════════════════════════
# 5. Visualization — Manhattan and QQ plots
# ═══════════════════════════════════════════════════════════════════════════

def manhattan_plot(trait: str, gwas_dir: str, output_dir: str,
                   bonferroni_p: float, suggestive_p: float = 1e-5):
    """Generate Manhattan plot for a single trait."""
    pval_file = os.path.join(gwas_dir, f"gwas_pvals_{trait}.csv.gz")
    if not os.path.exists(pval_file):
        return

    df = pd.read_csv(pval_file)
    df = df.dropna(subset=["P"])
    df = df[df["P"] > 0]
    df["-logP"] = -np.log10(df["P"])

    # Assign chromosome order
    chrom_str = df["CHROM"].astype(str)
    # Extract numeric part
    chrom_num = chrom_str.str.extract(r'(\d+)')[0].astype(float)
    df["CHROM_NUM"] = chrom_num

    # Remove non-standard chromosomes
    df = df.dropna(subset=["CHROM_NUM"])
    df["CHROM_NUM"] = df["CHROM_NUM"].astype(int)
    df = df.sort_values(["CHROM_NUM", "POS"])

    # Create cumulative position for x-axis
    chrom_groups = df.groupby("CHROM_NUM")
    chrom_offsets = {}
    cumulative = 0
    chrom_centers = {}

    for chrom in sorted(df["CHROM_NUM"].unique()):
        chrom_offsets[chrom] = cumulative
        sub = chrom_groups.get_group(chrom)
        chrom_centers[chrom] = cumulative + sub["POS"].median()
        cumulative += sub["POS"].max() + 1000000  # gap between chromosomes

    df["CUM_POS"] = df.apply(
        lambda r: r["POS"] + chrom_offsets.get(r["CHROM_NUM"], 0), axis=1
    )

    # Plot
    fig, ax = plt.subplots(figsize=(14, 5))
    colors = ["#2171B5", "#6BAED6"]

    for chrom in sorted(df["CHROM_NUM"].unique()):
        sub = df[df["CHROM_NUM"] == chrom]
        color = colors[int(chrom) % 2]
        ax.scatter(sub["CUM_POS"], sub["-logP"], s=2, color=color,
                   alpha=0.5, rasterized=True)

    # Threshold lines
    ax.axhline(-np.log10(bonferroni_p), color="red", linewidth=0.8,
               linestyle="--", label=f"Bonferroni (p={bonferroni_p:.1e})")
    ax.axhline(-np.log10(suggestive_p), color="blue", linewidth=0.8,
               linestyle=":", label=f"Suggestive (p={suggestive_p:.0e})")

    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("-log₁₀(p)", fontsize=12)
    ax.set_title(f"Manhattan Plot: {trait}", fontsize=13)
    ax.set_xticks(list(chrom_centers.values()))
    ax.set_xticklabels([str(c) for c in sorted(chrom_centers.keys())])
    ax.legend(fontsize=8)
    ax.set_xlim(0, cumulative)
    ax.set_ylim(0, max(df["-logP"].max() + 1, -np.log10(bonferroni_p) + 1))

    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f"manhattan_{trait}.png"), dpi=150)
    plt.close()


def qq_plot(trait: str, gwas_dir: str, output_dir: str):
    """Generate QQ plot for a single trait."""
    pval_file = os.path.join(gwas_dir, f"gwas_pvals_{trait}.csv.gz")
    if not os.path.exists(pval_file):
        return

    df = pd.read_csv(pval_file)
    pvals = df["P"].dropna().values
    pvals = pvals[pvals > 0]
    pvals = np.sort(pvals)

    n = len(pvals)
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    observed = -np.log10(pvals)

    # Compute genomic inflation factor (lambda)
    chi2 = stats.chi2.isf(pvals, df=1)
    lambda_gc = np.median(chi2) / stats.chi2.ppf(0.5, df=1)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(expected, observed, s=2, color="steelblue", alpha=0.5, rasterized=True)
    max_val = max(expected.max(), observed.max()) + 0.5
    ax.plot([0, max_val], [0, max_val], "r--", linewidth=1)
    ax.set_xlabel("Expected -log₁₀(p)", fontsize=12)
    ax.set_ylabel("Observed -log₁₀(p)", fontsize=12)
    ax.set_title(f"QQ Plot: {trait} (λ = {lambda_gc:.3f})", fontsize=13)
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)

    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f"qq_{trait}.png"), dpi=150)
    plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# 6. Main runner
# ═══════════════════════════════════════════════════════════════════════════

class GWASMLHybrid:
    """Run GWAS, annotate genes, and cross-reference with ML predictions."""

    def __init__(
        self,
        phenotype_file: str,
        vcf_file: str,
        ml_results_file: str,
        gff_file: Optional[str] = None,
        accession_col: str = "accession_id",
        exclude_cols: Optional[list] = None,
        n_pop_pcs: int = 5,
        ml_threshold: float = 0.3,
        output_dir: str = "gwas_results",
        remove_outliers: bool = True,
    ):
        self.phenotype_file = phenotype_file
        self.vcf_file = vcf_file
        self.ml_results_file = ml_results_file
        self.gff_file = gff_file
        self.accession_col = accession_col
        self.exclude_cols = exclude_cols or ["replicate_id"]
        self.n_pop_pcs = n_pop_pcs
        self.ml_threshold = ml_threshold
        self.output_dir = output_dir
        self.remove_outliers = remove_outliers

    def run(self):
        os.makedirs(self.output_dir, exist_ok=True)

        # ── Load data ─────────────────────────────────────────────────
        pheno_df, traits = load_phenotype(
            self.phenotype_file, self.accession_col, self.exclude_cols,
            remove_outliers=self.remove_outliers
        )
        geno_matrix, matched_ids, snp_info = load_genotype(
            self.vcf_file, pheno_df[self.accession_col].tolist()
        )

        n_snps = geno_matrix.shape[1]
        bonferroni_p = 0.05 / n_snps

        # ── Population structure PCs ──────────────────────────────────
        pop_pcs = compute_population_pcs(geno_matrix, self.n_pop_pcs)

        # ── Run GWAS ──────────────────────────────────────────────────
        gwas_dir = os.path.join(self.output_dir, "per_trait")
        gwas_summary = run_gwas_all_traits(
            pheno_df, traits, geno_matrix, matched_ids,
            snp_info, pop_pcs, gwas_dir,
        )

        # ── Gene annotations ──────────────────────────────────────────
        if self.gff_file is None:
            self.gff_file = download_tair10_gff(
                os.path.join(self.output_dir, "TAIR10_GFF3_genes.gff")
            )
        gene_df = load_gene_annotations(self.gff_file)

        # ── Cross-reference with ML ───────────────────────────────────
        cross_dir = os.path.join(self.output_dir, "ml_cross_reference")
        highlight = cross_reference_gwas_ml(
            gwas_summary, self.ml_results_file, gwas_dir,
            gene_df, cross_dir, self.ml_threshold,
        )

        # ── Plots (for top ML traits) ────────────────────────────────
        ml_df = pd.read_csv(self.ml_results_file)
        top_traits = ml_df.nlargest(20, "Best_r")["Trait"].tolist()

        plot_dir = os.path.join(self.output_dir, "plots")
        logging.info(f"Generating Manhattan and QQ plots for top {len(top_traits)} traits")
        for trait in top_traits:
            manhattan_plot(trait, gwas_dir, plot_dir, bonferroni_p)
            qq_plot(trait, gwas_dir, plot_dir)

        # ── Print summary ─────────────────────────────────────────────
        self._print_summary(gwas_summary, ml_df, highlight, bonferroni_p)

        logging.info(f"All GWAS-ML hybrid results saved to {self.output_dir}/")

    def _print_summary(self, gwas_summary, ml_df, highlight, bonferroni_p):
        merged = gwas_summary.merge(ml_df[["Trait", "Best_r", "Best_Model"]],
                                     on="Trait", how="left")

        print("\n" + "=" * 75)
        print("GWAS-ML HYBRID RESULTS SUMMARY")
        print("=" * 75)
        print(f"  Bonferroni threshold: p < {bonferroni_p:.2e}")
        print(f"  Suggestive threshold: p < 1e-05")
        print(f"  Total traits analyzed: {len(gwas_summary)}")
        print(f"  Traits with Bonferroni hits: "
              f"{(gwas_summary['N_Bonferroni'] > 0).sum()}")
        print(f"  Traits with suggestive hits: "
              f"{(gwas_summary['N_Suggestive'] > 0).sum()}")

        print(f"\n  Top 15 traits (by ML performance + GWAS hits):")
        top = merged.nlargest(15, "Best_r")
        print(f"  {'Trait':<22s} {'ML_r':>6s} {'Model':<8s} "
              f"{'Bonf':>5s} {'Sugg':>5s} {'Min_P':>10s} {'Top_SNP':>12s}")
        print(f"  {'-'*72}")
        for _, r in top.iterrows():
            print(f"  {r['Trait']:<22s} {r['Best_r']:>6.3f} {r['Best_Model']:<8s} "
                  f"{int(r['N_Bonferroni']):>5d} {int(r['N_Suggestive']):>5d} "
                  f"{r['Min_P']:>10.2e} {str(r.get('Top_SNP','')):>12s}")

        if highlight is not None and not highlight.empty:
            print(f"\n  Highlight: Top gene annotations for ML-significant traits:")
            shown = set()
            for _, r in highlight.iterrows():
                key = f"{r['Trait']}_{r.get('Nearest_Gene_Name','')}"
                if key in shown:
                    continue
                shown.add(key)
                if len(shown) > 20:
                    break
                print(f"    {r['Trait']:<20s} p={r['P_value']:.2e}  "
                      f"{r.get('Nearest_Gene_Name','N/A'):<15s} "
                      f"{str(r.get('Gene_Description',''))[:50]}")

        print("=" * 75)


# ═══════════════════════════════════════════════════════════════════════════
# CLI entry point
# ═══════════════════════════════════════════════════════════════════════════
def parse_args():
    """Parse command-line arguments."""
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "GWAS-ML Hybrid Pipeline: per-SNP association testing with\n"
            "population structure correction, cross-referenced with ML\n"
            "prediction performance, and automatic gene annotation."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python gwas_ml_hybrid.py \\\n"
            "      --phenotype data/phenotypes.csv \\\n"
            "      --vcf data/genotypes.vcf \\\n"
            "      --ml-results results/phenotype_prediction_results.csv\n\n"
            "  python gwas_ml_hybrid.py \\\n"
            "      --phenotype data/rice_traits.csv \\\n"
            "      --vcf data/rice_44k.vcf \\\n"
            "      --ml-results results_rice/phenotype_prediction_results.csv \\\n"
            "      --gff data/rice_MSU7_genes.gff3 \\\n"
            "      --accession-col HybID \\\n"
            "      --exclude-cols replicate_id batch \\\n"
            "      --n-pop-pcs 5 \\\n"
            "      --ml-threshold 0.3 \\\n"
            "      --output gwas_rice/\n"
        ),
    )

    # ── Required arguments ────────────────────────────────────────────
    parser.add_argument(
        "--phenotype", required=True,
        help="Path to phenotype CSV file (rows=accessions, columns=traits)",
    )
    parser.add_argument(
        "--vcf", required=True,
        help="Path to VCF genotype file (biallelic SNPs)",
    )
    parser.add_argument(
        "--ml-results", required=True,
        help=(
            "Path to ML results CSV from adaptive_ensemble_pipeline.py "
            "(phenotype_prediction_results.csv)"
        ),
    )

    # ── Optional arguments ────────────────────────────────────────────
    parser.add_argument(
        "--gff", default=None,
        help=(
            "Path to GFF3 gene annotation file for nearest-gene annotation. "
            "If not provided, GWAS runs without gene annotation."
        ),
    )
    parser.add_argument(
        "--output", default="gwas_results",
        help="Output directory for all GWAS results (default: gwas_results/)",
    )
    parser.add_argument(
        "--accession-col", default="accession_id",
        help="Column name for accession IDs in phenotype CSV (default: accession_id)",
    )
    parser.add_argument(
        "--exclude-cols", nargs="*", default=["replicate_id"],
        help="Non-trait columns to exclude (default: replicate_id)",
    )
    parser.add_argument(
        "--n-pop-pcs", type=int, default=5,
        help="Number of PCs for population structure correction (default: 5)",
    )
    parser.add_argument(
        "--ml-threshold", type=float, default=0.3,
        help=(
            "Minimum ML Pearson r to include a trait in GWAS cross-reference "
            "(default: 0.3)"
        ),
    )
    parser.add_argument(
        "--no-outlier-removal", action="store_true", default=False,
        help=(
            "Disable per-trait outlier removal (|z|>3). Recommended when your "
            "phenotype data is already cleaned or when sample size is small "
            "(n < 200) and you cannot afford to lose accessions."
        ),
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    pipeline = GWASMLHybrid(
        phenotype_file=args.phenotype,
        vcf_file=args.vcf,
        ml_results_file=args.ml_results,
        gff_file=args.gff,
        accession_col=args.accession_col,
        exclude_cols=args.exclude_cols,
        n_pop_pcs=args.n_pop_pcs,
        ml_threshold=args.ml_threshold,
        output_dir=args.output,
        remove_outliers=not args.no_outlier_removal,
    )
    pipeline.run()
