"""
Adaptive Ensemble ML Pipeline for Genotype-to-Phenotype Prediction.

Automatically scales model complexity, PCA components, CV strategy, and
ensemble composition based on dataset size. Works for datasets from
~30 to 100,000+ samples.

Usage (command line):
    python adaptive_ensemble_pipeline.py \\
        --phenotype data/phenotypes.csv \\
        --vcf data/genotypes.vcf \\
        --output results/

    python adaptive_ensemble_pipeline.py \\
        --phenotype data/rice_traits.csv \\
        --vcf data/rice_44k.vcf \\
        --accession-col HybID \\
        --exclude-cols replicate_id batch \\
        --output results_rice/

    For full options:
        python adaptive_ensemble_pipeline.py --help

Usage (Python API):
    from adaptive_ensemble_pipeline import GenomePhenotypePredictor

    predictor = GenomePhenotypePredictor(
        phenotype_file="pheno.csv",
        vcf_file="geno.vcf",
        accession_col="accession_id",
        exclude_cols=["replicate_id"],
    )
    predictor.run()
"""

import pandas as pd
import numpy as np
import allel
import json
from sklearn.ensemble import RandomForestRegressor, StackingRegressor
from sklearn.linear_model import Ridge, ElasticNet
from sklearn.multioutput import MultiOutputRegressor
from sklearn.model_selection import (
    train_test_split,
    cross_val_score,
    KFold,
    RepeatedKFold,
)
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.impute import KNNImputer
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import Pipeline
from sklearn.base import clone, BaseEstimator, RegressorMixin
from scipy.stats import pearsonr
from sklearn.metrics import r2_score, mean_squared_error
from xgboost import XGBRegressor
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import logging
import os
import sys
import warnings
from dataclasses import dataclass, field
from typing import Optional

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)


# ═══════════════════════════════════════════════════════════════════════════
# Adaptive configuration — everything scales from sample count
# ═══════════════════════════════════════════════════════════════════════════
@dataclass
class AdaptiveConfig:
    """Automatically determines all hyperparameters from dataset dimensions."""

    n_samples: int
    n_features: int
    n_traits: int

    # Computed fields (populated by __post_init__)
    regime: str = field(init=False)
    n_pca: int = field(init=False)
    test_size: float = field(init=False)
    cv_folds: int = field(init=False)
    cv_repeats: int = field(init=False)
    use_nn: bool = field(init=False)
    use_stacking: bool = field(init=False)
    rf_params: dict = field(init=False)
    xgb_params: dict = field(init=False)
    nn_layers: tuple = field(init=False)
    nn_max_iter: int = field(init=False)
    ridge_alpha: float = field(init=False)
    stacking_cv: int = field(init=False)
    knn_neighbors: int = field(init=False)
    variance_threshold: float = field(init=False)
    min_cv_train: int = field(init=False)

    def __post_init__(self):
        n = self.n_samples

        # ── Regime classification ──────────────────────────────────────
        if n < 50:
            self.regime = "tiny"
        elif n < 150:
            self.regime = "small"
        elif n < 500:
            self.regime = "medium"
        elif n < 2000:
            self.regime = "large"
        else:
            self.regime = "xlarge"

        # ── Test split size ────────────────────────────────────────────
        # Larger test fraction for tiny data so test set has ≥8 samples
        if n < 50:
            self.test_size = 0.25
        elif n < 200:
            self.test_size = 0.20
        else:
            self.test_size = 0.15

        n_train = int(n * (1 - self.test_size))

        # ── Cross-validation ───────────────────────────────────────────
        if n < 50:
            self.cv_folds = 3
            self.cv_repeats = 5       # repeat more to stabilize with tiny data
        elif n < 150:
            self.cv_folds = 5
            self.cv_repeats = 3
        elif n < 500:
            self.cv_folds = 10
            self.cv_repeats = 1
        else:
            self.cv_folds = 10
            self.cv_repeats = 1

        # Minimum training samples in any CV fold (used to cap PCA)
        self.min_cv_train = int(n_train * (1 - 1 / self.cv_folds))

        # ── PCA components ─────────────────────────────────────────────
        # Rule: ~70-80% of the min CV training fold, capped at 500
        pca_ceiling = max(10, int(self.min_cv_train * 0.75))
        self.n_pca = min(pca_ceiling, self.n_features, 500)

        # ── Feature pre-filtering ──────────────────────────────────────
        self.variance_threshold = 0.01 if self.n_features > 10000 else 0.0

        # ── KNN imputer ───────────────────────────────────────────────
        self.knn_neighbors = min(5, max(2, n // 20))

        # ── Model selection ────────────────────────────────────────────
        self.use_nn = (n >= 60)       # NN needs enough data
        self.use_stacking = (n >= 80) # Stacking needs internal CV headroom

        # ── Random Forest ──────────────────────────────────────────────
        if n < 100:
            self.rf_params = dict(
                n_estimators=100, max_depth=3, min_samples_leaf=5,
                max_features="sqrt", random_state=RANDOM_STATE,
            )
        elif n < 500:
            self.rf_params = dict(
                n_estimators=200, max_depth=5, min_samples_leaf=3,
                max_features="sqrt", random_state=RANDOM_STATE,
            )
        else:
            self.rf_params = dict(
                n_estimators=500, max_depth=8, min_samples_leaf=2,
                max_features="sqrt", random_state=RANDOM_STATE,
            )

        # ── XGBoost ───────────────────────────────────────────────────
        if n < 100:
            self.xgb_params = dict(
                n_estimators=100, max_depth=2, learning_rate=0.05,
                subsample=0.8, colsample_bytree=0.8, reg_alpha=1.0,
                reg_lambda=2.0, random_state=RANDOM_STATE, verbosity=0,
            )
        elif n < 500:
            self.xgb_params = dict(
                n_estimators=200, max_depth=4, learning_rate=0.05,
                subsample=0.8, colsample_bytree=0.8, reg_alpha=0.5,
                reg_lambda=1.5, random_state=RANDOM_STATE, verbosity=0,
            )
        else:
            self.xgb_params = dict(
                n_estimators=500, max_depth=6, learning_rate=0.05,
                subsample=0.8, colsample_bytree=0.8, reg_alpha=0.1,
                reg_lambda=1.0, random_state=RANDOM_STATE, verbosity=0,
            )

        # ── Neural Network ─────────────────────────────────────────────
        if n < 100:
            self.nn_layers = (32,)
            self.nn_max_iter = 300
        elif n < 500:
            self.nn_layers = (64, 32)
            self.nn_max_iter = 500
        elif n < 2000:
            self.nn_layers = (128, 64)
            self.nn_max_iter = 500
        else:
            self.nn_layers = (256, 128, 64)
            self.nn_max_iter = 1000

        # ── Ridge (simple baseline / stacking meta-learner) ───────────
        self.ridge_alpha = 10.0 if n < 200 else 1.0

        # ── Stacking internal CV ──────────────────────────────────────
        self.stacking_cv = 3 if n < 200 else 5

    def summary(self) -> str:
        return (
            f"AdaptiveConfig(regime={self.regime}, n={self.n_samples}, "
            f"features={self.n_features}, pca={self.n_pca}, "
            f"cv={self.cv_folds}×{self.cv_repeats}, test={self.test_size}, "
            f"nn={self.use_nn}, stack={self.use_stacking})"
        )


# ═══════════════════════════════════════════════════════════════════════════
# Single-output MLP wrapper (sklearn ≥1.6 compatible)
# ═══════════════════════════════════════════════════════════════════════════
class SingleOutputMLP(BaseEstimator, RegressorMixin):
    """Wraps MLPRegressor to guarantee 1-D output for StackingRegressor."""

    _estimator_type = "regressor"

    def __init__(self, hidden_layer_sizes=(64, 32), max_iter=500,
                 random_state=42, early_stopping=True,
                 validation_fraction=0.15, n_iter_no_change=20,
                 alpha=0.001):
        self.hidden_layer_sizes = hidden_layer_sizes
        self.max_iter = max_iter
        self.random_state = random_state
        self.early_stopping = early_stopping
        self.validation_fraction = validation_fraction
        self.n_iter_no_change = n_iter_no_change
        self.alpha = alpha

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.estimator_type = "regressor"
        try:
            from sklearn.utils._tags import RegressorTags
            tags.regressor_tags = RegressorTags()
        except ImportError:
            pass
        tags.target_tags.required = True
        return tags

    def fit(self, X, y):
        self.model_ = MLPRegressor(
            hidden_layer_sizes=self.hidden_layer_sizes,
            max_iter=self.max_iter,
            random_state=self.random_state,
            early_stopping=self.early_stopping,
            validation_fraction=self.validation_fraction,
            n_iter_no_change=self.n_iter_no_change,
            alpha=self.alpha,
        )
        self.model_.fit(X, y)
        return self

    def predict(self, X):
        return self.model_.predict(X).ravel()


# ═══════════════════════════════════════════════════════════════════════════
# Main pipeline
# ═══════════════════════════════════════════════════════════════════════════
class GenomePhenotypePredictor:
    """Adaptive genotype→phenotype prediction pipeline.

    Automatically configures model complexity, PCA dimensionality,
    cross-validation strategy, and ensemble composition based on the
    actual dataset dimensions.
    """

    def __init__(
        self,
        phenotype_file: str,
        vcf_file: str,
        accession_col: str = "accession_id",
        exclude_cols: Optional[list] = None,
        n_components: Optional[int] = None,   # None = auto
        test_size: Optional[float] = None,     # None = auto
        output_dir: str = "results",
        remove_outliers: bool = True,
    ):
        self.phenotype_file = phenotype_file
        self.vcf_file = vcf_file
        self.accession_col = accession_col
        self.exclude_cols = exclude_cols or ["replicate_id"]
        self._user_n_components = n_components
        self._user_test_size = test_size
        self.output_dir = output_dir
        self.remove_outliers = remove_outliers

        # Populated during run()
        self.phenotype_df = None
        self.genotype_df = None
        self.target_traits = None
        self.config = None

    # ──────────────────────────────────────────────────────────────────
    # 1. Phenotype loading & cleaning
    # ──────────────────────────────────────────────────────────────────
    def load_phenotype_data(self) -> tuple:
        logging.info(f"Loading phenotype data from {self.phenotype_file}")
        df = pd.read_csv(self.phenotype_file)

        if self.accession_col not in df.columns:
            raise ValueError(
                f"Accession column '{self.accession_col}' not found. "
                f"Available columns: {list(df.columns[:20])}"
            )

        # Identify trait columns (everything except accession + exclude)
        skip = set([self.accession_col] + self.exclude_cols)
        self.target_traits = [c for c in df.columns if c not in skip]
        logging.info(f"Detected {len(self.target_traits)} trait columns")

        # Aggregate replicates
        df_agg = df.groupby(self.accession_col)[self.target_traits].mean().reset_index()
        logging.info(f"After aggregation: {len(df_agg)} unique accessions")

        # KNN imputation (adaptive neighbors)
        k = min(5, max(2, len(df_agg) // 20))
        imputer = KNNImputer(n_neighbors=k)
        df_agg[self.target_traits] = imputer.fit_transform(df_agg[self.target_traits])

        # Robust scaling
        scaler = RobustScaler()
        df_agg[self.target_traits] = scaler.fit_transform(df_agg[self.target_traits])

        # Outlier removal (|z| > 3), trait-specific
        if self.remove_outliers:
            before = len(df_agg)
            for trait in self.target_traits:
                mu, sigma = df_agg[trait].mean(), df_agg[trait].std()
                if sigma > 0:
                    df_agg = df_agg[np.abs((df_agg[trait] - mu) / sigma) <= 3]
            logging.info(f"Outlier removal: {before} → {len(df_agg)} accessions")
        else:
            logging.info(f"Outlier removal: SKIPPED (--no-outlier-removal), keeping {len(df_agg)} accessions")

        self.phenotype_df = df_agg
        return df_agg, self.target_traits

    # ──────────────────────────────────────────────────────────────────
    # 2. Genotype loading & robust accession matching
    # ──────────────────────────────────────────────────────────────────
    def load_genotype_data(self, phenotype_accession_ids: list) -> pd.DataFrame:
        logging.info(f"Loading genotype data from {self.vcf_file}")
        callset = allel.read_vcf(
            self.vcf_file,
            fields=["samples", "calldata/GT", "variants/ID"],
        )
        if callset is None:
            raise ValueError("VCF file could not be read")

        genotypes = callset["calldata/GT"].sum(axis=2).astype(float)
        sample_ids = list(callset["samples"])
        snp_ids = list(callset["variants/ID"])
        logging.info(f"VCF: {len(sample_ids)} samples, {len(snp_ids)} SNPs")

        genotypes = genotypes.T
        genotypes[genotypes < 0] = np.nan

        if len(snp_ids) != len(set(snp_ids)):
            logging.warning("Duplicate SNP IDs — generating unique IDs")
            snp_ids = [f"SNP_{i}" for i in range(len(snp_ids))]

        genotype_df = pd.DataFrame(
            genotypes,
            index=[str(s) for s in sample_ids],
            columns=[str(s) for s in snp_ids],
        )

        # ── Multi-strategy accession matching ──────────────────────────
        pheno_ids = set(str(a) for a in phenotype_accession_ids)
        vcf_ids = set(genotype_df.index)

        # Strategy A: direct string match
        direct = pheno_ids & vcf_ids
        logging.info(f"ID matching — direct: {len(direct)} matches")

        if len(direct) >= max(10, len(pheno_ids) * 0.3):
            genotype_df = genotype_df.loc[genotype_df.index.isin(pheno_ids)]
        else:
            # Strategy B: numeric match
            vcf_numeric, pheno_numeric = {}, {}
            for vid in vcf_ids:
                try:
                    vcf_numeric[str(int(vid))] = vid
                except ValueError:
                    continue
            for pid in pheno_ids:
                try:
                    pheno_numeric[str(int(pid))] = pid
                except ValueError:
                    continue

            numeric_matches = set(vcf_numeric) & set(pheno_numeric)
            logging.info(f"ID matching — numeric: {len(numeric_matches)} matches")

            if len(numeric_matches) >= max(10, len(pheno_ids) * 0.3):
                keep = [vcf_numeric[k] for k in numeric_matches]
                genotype_df = genotype_df.loc[keep]
                rename_map = {vcf_numeric[k]: pheno_numeric[k] for k in numeric_matches}
                genotype_df.index = genotype_df.index.map(rename_map)
            else:
                # Strategy C: common prefix stripping
                sample_vcf = list(vcf_ids)[:50]
                prefix_len = 0
                for plen in range(1, 6):
                    prefixes = set(s[:plen] for s in sample_vcf if len(s) > plen)
                    if len(prefixes) == 1:
                        prefix_len = plen
                    else:
                        break

                if prefix_len > 0:
                    prefix = list(set(s[:prefix_len] for s in sample_vcf))[0]
                    logging.info(f"ID matching — trying prefix strip: '{prefix}'")
                    vcf_stripped = {}
                    for vid in vcf_ids:
                        if vid.startswith(prefix):
                            try:
                                vcf_stripped[str(int(vid[prefix_len:]))] = vid
                            except ValueError:
                                continue

                    stripped_matches = set(vcf_stripped) & set(pheno_numeric)
                    logging.info(f"ID matching — prefix-stripped: {len(stripped_matches)}")

                    if len(stripped_matches) >= 10:
                        keep = [vcf_stripped[k] for k in stripped_matches]
                        genotype_df = genotype_df.loc[keep]
                        rename_map = {vcf_stripped[k]: pheno_numeric[k]
                                      for k in stripped_matches}
                        genotype_df.index = genotype_df.index.map(rename_map)
                    else:
                        raise ValueError(
                            f"Could not match accession IDs.\n"
                            f"  VCF IDs (first 10): {sample_ids[:10]}\n"
                            f"  Phenotype IDs (first 10): {list(pheno_ids)[:10]}\n"
                            f"  Direct: {len(direct)}, Numeric: {len(numeric_matches)}, "
                            f"Stripped: {len(stripped_matches)}"
                        )
                else:
                    raise ValueError(
                        f"Could not match IDs.\n"
                        f"  VCF (first 10): {sample_ids[:10]}\n"
                        f"  Pheno (first 10): {list(pheno_ids)[:10]}"
                    )

        if genotype_df.empty:
            raise ValueError("No overlapping accessions.")

        # Impute missing genotypes
        k = min(5, max(2, len(genotype_df) // 20))
        imp = KNNImputer(n_neighbors=k)
        genotype_df = pd.DataFrame(
            imp.fit_transform(genotype_df),
            index=genotype_df.index,
            columns=genotype_df.columns,
        )

        self.genotype_df = genotype_df
        logging.info(f"Genotype data: {genotype_df.shape[0]} accessions, "
                     f"{genotype_df.shape[1]} SNPs")
        return genotype_df

    # ──────────────────────────────────────────────────────────────────
    # 3. Merge & pre-filter
    # ──────────────────────────────────────────────────────────────────
    def merge_data(self) -> tuple:
        logging.info("Merging phenotype and genotype data")
        pheno = self.phenotype_df.copy()
        pheno[self.accession_col] = pheno[self.accession_col].astype(str)
        self.genotype_df.index = self.genotype_df.index.astype(str)

        merged = pheno.merge(
            self.genotype_df,
            left_on=self.accession_col,
            right_index=True,
            how="inner",
        )
        if merged.empty:
            raise ValueError("No overlapping accessions after merge.")

        logging.info(f"Merged: {len(merged)} accessions")

        feature_cols = [c for c in merged.columns
                        if c not in self.target_traits and c != self.accession_col]
        X = merged[feature_cols].astype(float)
        y = merged[self.target_traits].astype(float)

        # Variance-based pre-filtering for high-dimensional genotype data
        if X.shape[1] > 10000:
            before = X.shape[1]
            vt = VarianceThreshold(threshold=0.01)
            X_filtered = vt.fit_transform(X)
            kept_cols = X.columns[vt.get_support()]
            X = pd.DataFrame(X_filtered, index=X.index, columns=kept_cols)
            logging.info(f"Variance filter: {before} → {X.shape[1]} features")

        logging.info(f"Features: {X.shape[1]}  |  Targets: {y.shape[1]}")

        # ── Build adaptive config ─────────────────────────────────────
        self.config = AdaptiveConfig(
            n_samples=X.shape[0],
            n_features=X.shape[1],
            n_traits=y.shape[1],
        )
        # Allow user overrides
        if self._user_n_components is not None:
            self.config.n_pca = min(self._user_n_components, self.config.min_cv_train)
        if self._user_test_size is not None:
            self.config.test_size = self._user_test_size

        logging.info(f"Config: {self.config.summary()}")
        return X, y

    # ──────────────────────────────────────────────────────────────────
    # 4. Pipeline builder
    # ──────────────────────────────────────────────────────────────────
    def _make_pipeline(self, estimator) -> Pipeline:
        return Pipeline([
            ("scaler", StandardScaler()),
            ("pca", PCA(n_components=self.config.n_pca)),
            ("model", estimator),
        ])

    # ──────────────────────────────────────────────────────────────────
    # 5. Adaptive model builder
    # ──────────────────────────────────────────────────────────────────
    def _build_models(self) -> dict:
        """Build model dict based on adaptive config."""
        cfg = self.config
        models = {}

        # Always include Ridge as a simple baseline
        models["Ridge"] = self._make_pipeline(
            MultiOutputRegressor(Ridge(alpha=cfg.ridge_alpha))
        )

        # Always include RF
        models["RF"] = self._make_pipeline(
            MultiOutputRegressor(RandomForestRegressor(**cfg.rf_params))
        )

        # Always include XGBoost
        models["XGB"] = self._make_pipeline(
            MultiOutputRegressor(XGBRegressor(**cfg.xgb_params))
        )

        # NN only if enough data
        if cfg.use_nn:
            models["NN"] = self._make_pipeline(
                MultiOutputRegressor(
                    SingleOutputMLP(
                        hidden_layer_sizes=cfg.nn_layers,
                        max_iter=cfg.nn_max_iter,
                        random_state=RANDOM_STATE,
                        alpha=0.01 if cfg.n_samples < 200 else 0.001,
                    )
                )
            )
            logging.info(f"NN included: layers={cfg.nn_layers}")
        else:
            logging.info(f"NN skipped: too few samples ({cfg.n_samples})")

        return models

    # ──────────────────────────────────────────────────────────────────
    # 6. Training
    # ──────────────────────────────────────────────────────────────────
    def train_ensemble(self, X: pd.DataFrame, y: pd.DataFrame) -> dict:
        logging.info("Training ensemble models")
        cfg = self.config

        # Train/test split
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=cfg.test_size, random_state=RANDOM_STATE,
        )
        logging.info(f"Train: {len(X_train)}  |  Test: {len(X_test)}")

        # Build & fit base models
        models = self._build_models()
        predictions = {}

        for name, pipe in models.items():
            logging.info(f"  Training {name}...")
            pipe.fit(X_train, y_train)
            predictions[name] = pipe.predict(X_test)
            logging.info(f"  {name} done")

        # ── Stacking (per-trait) ──────────────────────────────────────
        if cfg.use_stacking:
            logging.info("  Training Stacking...")
            scaler_pca = Pipeline([
                ("scaler", StandardScaler()),
                ("pca", PCA(n_components=cfg.n_pca)),
            ])
            X_train_pca = scaler_pca.fit_transform(X_train)
            X_test_pca = scaler_pca.transform(X_test)

            stack_preds = np.zeros((len(X_test), len(self.target_traits)))
            stack_models = []

            # Build base estimator list for stacking
            base_estimators = [
                ("rf", RandomForestRegressor(**cfg.rf_params)),
                ("xgb", XGBRegressor(**cfg.xgb_params)),
            ]
            if cfg.use_nn:
                base_estimators.append(
                    ("nn", SingleOutputMLP(
                        hidden_layer_sizes=cfg.nn_layers,
                        max_iter=cfg.nn_max_iter,
                        random_state=RANDOM_STATE,
                        alpha=0.01 if cfg.n_samples < 200 else 0.001,
                    ))
                )

            # Use Ridge as meta-learner (more stable than RF for small data)
            meta_learner = Ridge(alpha=cfg.ridge_alpha)

            for i, trait in enumerate(self.target_traits):
                stack = StackingRegressor(
                    estimators=[(n, clone(e)) for n, e in base_estimators],
                    final_estimator=clone(meta_learner),
                    cv=cfg.stacking_cv,
                )
                stack.fit(X_train_pca, y_train.iloc[:, i])
                stack_models.append(stack)
                stack_preds[:, i] = stack.predict(X_test_pca)

            models["Stack"] = stack_models
            predictions["Stack"] = stack_preds
            logging.info("  Stacking done")
        else:
            logging.info(f"  Stacking skipped: too few samples ({cfg.n_samples})")

        # ── Cross-validation ──────────────────────────────────────────
        logging.info("Running cross-validation...")
        if cfg.cv_repeats > 1:
            cv = RepeatedKFold(
                n_splits=cfg.cv_folds, n_repeats=cfg.cv_repeats,
                random_state=RANDOM_STATE,
            )
        else:
            cv = KFold(
                n_splits=cfg.cv_folds, shuffle=True,
                random_state=RANDOM_STATE,
            )

        cv_scores = {}
        for name, pipe in models.items():
            if name == "Stack":
                # Per-trait stacking CV
                trait_scores = []
                for i, trait in enumerate(self.target_traits):
                    stack_pipe = self._make_pipeline(
                        StackingRegressor(
                            estimators=[(n, clone(e)) for n, e in base_estimators],
                            final_estimator=clone(meta_learner),
                            cv=cfg.stacking_cv,
                        )
                    )
                    score = cross_val_score(
                        stack_pipe, X, y.iloc[:, i], cv=cv, scoring="r2"
                    ).mean()
                    trait_scores.append(score)
                cv_scores[name] = trait_scores
                logging.info(f"  CV R² {name}: {np.mean(trait_scores):.4f} (mean)")
            else:
                scores = cross_val_score(pipe, X, y, cv=cv, scoring="r2")
                cv_scores[name] = scores.mean()
                logging.info(f"  CV R² {name}: {scores.mean():.4f}")

        return {
            "models": models,
            "predictions": predictions,
            "cv_scores": cv_scores,
            "splits": {
                "X_train": X_train, "X_test": X_test,
                "y_train": y_train, "y_test": y_test,
            },
        }

    # ──────────────────────────────────────────────────────────────────
    # 7. Evaluation
    # ──────────────────────────────────────────────────────────────────
    def evaluate(self, y_test, predictions: dict) -> dict:
        results = {}
        for model_name, preds in predictions.items():
            pearson = {}
            for i, trait in enumerate(self.target_traits):
                obs = y_test.iloc[:, i].values
                pred = preds[:, i]
                if np.std(obs) == 0 or np.std(pred) == 0:
                    pearson[trait] = 0.0
                else:
                    r, p = pearsonr(obs, pred)
                    pearson[trait] = r
            results[model_name] = pearson
        return results

    # ──────────────────────────────────────────────────────────────────
    # 8. RF uncertainty (from individual trees)
    # ──────────────────────────────────────────────────────────────────
    @staticmethod
    def uncertainty_from_rf(rf_pipe, X_test) -> tuple:
        X_t = X_test
        for name, step in rf_pipe.steps[:-1]:
            X_t = step.transform(X_t)
        multi_rf = rf_pipe.named_steps["model"]

        n_traits = len(multi_rf.estimators_)
        n_samples = X_t.shape[0]
        means = np.zeros((n_samples, n_traits))
        stds = np.zeros((n_samples, n_traits))

        for i, est in enumerate(multi_rf.estimators_):
            tree_preds = np.array([t.predict(X_t) for t in est.estimators_])
            means[:, i] = tree_preds.mean(axis=0)
            stds[:, i] = tree_preds.std(axis=0)
        return means, stds

    # ──────────────────────────────────────────────────────────────────
    # 9. Feature importance
    # ──────────────────────────────────────────────────────────────────
    def get_feature_importance(self, rf_pipe, top_n: int = 10) -> dict:
        pca = rf_pipe.named_steps["pca"]
        multi_rf = rf_pipe.named_steps["model"]
        pc_names = [f"PC_{j+1}" for j in range(pca.n_components_)]

        result = {}
        for i, trait in enumerate(self.target_traits):
            imp = multi_rf.estimators_[i].feature_importances_
            df = pd.DataFrame({"Feature": pc_names, "Importance": imp})
            result[trait] = df.sort_values("Importance", ascending=False).head(top_n)
        return result

    # ──────────────────────────────────────────────────────────────────
    # 10. Plotting
    # ──────────────────────────────────────────────────────────────────
    def plot_pearson_bars(self, pearson_dict: dict, model_name: str, save_dir: str):
        os.makedirs(save_dir, exist_ok=True)
        traits = list(pearson_dict.keys())
        vals = list(pearson_dict.values())
        order = np.argsort(vals)[::-1]
        traits = [traits[i] for i in order]
        vals = [vals[i] for i in order]
        colors = ["steelblue" if v >= 0 else "salmon" for v in vals]

        fig_width = max(10, len(traits) * 0.25)
        plt.figure(figsize=(fig_width, 6))
        plt.bar(range(len(traits)), vals, color=colors)
        plt.xticks(range(len(traits)), traits, rotation=90, fontsize=7)
        plt.xlabel("Trait")
        plt.ylabel("Pearson r")
        plt.title(f"Pearson Correlation by Trait ({model_name})")
        plt.axhline(0, color="black", linewidth=0.5)
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, f"pearson_correlations_{model_name}.png"),
                    dpi=150)
        plt.close()

    def plot_model_comparison(self, pearson_results: dict, save_dir: str):
        """Side-by-side bar plot comparing all models."""
        os.makedirs(save_dir, exist_ok=True)
        model_names = list(pearson_results.keys())
        traits = self.target_traits

        # Sort traits by best model performance
        best_r = [max(pearson_results[m].get(t, 0) for m in model_names) for t in traits]
        order = np.argsort(best_r)[::-1]
        traits_sorted = [traits[i] for i in order]

        x = np.arange(len(traits_sorted))
        width = 0.8 / len(model_names)
        colors = ["steelblue", "darkorange", "seagreen", "crimson", "mediumpurple"]

        fig_width = max(12, len(traits_sorted) * 0.3)
        plt.figure(figsize=(fig_width, 7))
        for j, mname in enumerate(model_names):
            vals = [pearson_results[mname].get(t, 0) for t in traits_sorted]
            plt.bar(x + j * width, vals, width, label=mname,
                    color=colors[j % len(colors)], alpha=0.8)

        plt.xticks(x + width * len(model_names) / 2, traits_sorted,
                   rotation=90, fontsize=6)
        plt.xlabel("Trait")
        plt.ylabel("Pearson r")
        plt.title("Model Comparison: Pearson r by Trait")
        plt.legend()
        plt.axhline(0, color="black", linewidth=0.5)
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, "model_comparison.png"), dpi=150)
        plt.close()
        logging.info("Saved model comparison plot")

    def plot_predicted_vs_observed(self, y_test, predictions: dict,
                                   pearson_results: dict, save_dir: str):
        os.makedirs(save_dir, exist_ok=True)

        # Determine best model per trait
        model_names = list(predictions.keys())
        colors = {"Ridge": "gray", "RF": "steelblue", "XGB": "darkorange",
                  "NN": "seagreen", "Stack": "crimson"}

        # Plot top positive traits (by best model)
        best_per_trait = {}
        for trait in self.target_traits:
            best_m = max(model_names,
                         key=lambda m: pearson_results[m].get(trait, -999))
            best_per_trait[trait] = (best_m, pearson_results[best_m].get(trait, 0))

        pos_traits = [(t, m, r) for t, (m, r) in best_per_trait.items() if r > 0]
        pos_traits.sort(key=lambda x: x[2], reverse=True)

        for trait, _, _ in pos_traits[:30]:  # Top 30
            i = self.target_traits.index(trait)
            obs = y_test.iloc[:, i].values

            fig, ax = plt.subplots(figsize=(7, 6))
            all_vals = list(obs)

            for mname in model_names:
                pred = predictions[mname][:, i]
                r = pearson_results[mname].get(trait, 0)
                c = colors.get(mname, "black")
                ax.scatter(obs, pred, label=f"{mname} (r={r:.3f})",
                           alpha=0.6, color=c, s=40)
                all_vals.extend(pred)

            vmin, vmax = min(all_vals), max(all_vals)
            margin = (vmax - vmin) * 0.05 + 1e-6
            ax.plot([vmin - margin, vmax + margin],
                    [vmin - margin, vmax + margin],
                    "k--", linewidth=1, label="Perfect")
            ax.set_xlabel("Observed")
            ax.set_ylabel("Predicted")
            ax.set_title(f"{trait}")
            ax.legend(fontsize=7)
            plt.tight_layout()
            plt.savefig(os.path.join(save_dir, f"{trait}.png"), dpi=150)
            plt.close()

        logging.info(f"Saved {min(len(pos_traits), 30)} scatter plots")

    # ──────────────────────────────────────────────────────────────────
    # 11. Save results
    # ──────────────────────────────────────────────────────────────────
    def save_results(self, y_test, predictions: dict, pearson_results: dict,
                     cv_scores: dict, uncertainty: dict):
        os.makedirs(self.output_dir, exist_ok=True)

        # ── predictions.csv ───────────────────────────────────────────
        frames = []
        for i, trait in enumerate(self.target_traits):
            cols = {f"{trait}_Observed": y_test.iloc[:, i].values}
            for mname, preds in predictions.items():
                cols[f"{trait}_{mname}_Predicted"] = preds[:, i]
                if mname in uncertainty:
                    _, stds = uncertainty[mname]
                    if stds.ndim == 2:
                        cols[f"{trait}_{mname}_Std"] = stds[:, i]
            frames.append(pd.DataFrame(cols))
        pred_df = pd.concat(frames, axis=1)
        pred_df.to_csv(os.path.join(self.output_dir, "predictions.csv"), index=False)

        # ── phenotype_prediction_results.csv ──────────────────────────
        rows = []
        for trait in self.target_traits:
            row = {"Trait": trait}
            for mname, pr in pearson_results.items():
                row[f"Pearson_r_{mname}"] = pr.get(trait, np.nan)
            # Best model for this trait
            best_m = max(pearson_results,
                         key=lambda m: pearson_results[m].get(trait, -999))
            row["Best_Model"] = best_m
            row["Best_r"] = pearson_results[best_m].get(trait, np.nan)
            rows.append(row)
        results_df = pd.DataFrame(rows)
        results_df.to_csv(
            os.path.join(self.output_dir, "phenotype_prediction_results.csv"),
            index=False,
        )

        # ── config.json ──────────────────────────────────────────────
        config_dict = {
            "regime": self.config.regime,
            "n_samples": self.config.n_samples,
            "n_features": self.config.n_features,
            "n_traits": self.config.n_traits,
            "n_pca": self.config.n_pca,
            "test_size": self.config.test_size,
            "cv_folds": self.config.cv_folds,
            "cv_repeats": self.config.cv_repeats,
            "use_nn": self.config.use_nn,
            "use_stacking": self.config.use_stacking,
        }
        with open(os.path.join(self.output_dir, "config.json"), "w") as f:
            json.dump(config_dict, f, indent=2)

        # ── Console summary ──────────────────────────────────────────
        self._print_summary(pearson_results, cv_scores)

        logging.info(f"All results saved to {self.output_dir}/")

    def _print_summary(self, pearson_results: dict, cv_scores: dict):
        cfg = self.config
        print("\n" + "=" * 70)
        print(f"RESULTS SUMMARY — Regime: {cfg.regime.upper()} "
              f"(n={cfg.n_samples}, PCA={cfg.n_pca}, "
              f"CV={cfg.cv_folds}×{cfg.cv_repeats})")
        print("=" * 70)

        # Per-model summary
        model_names = list(pearson_results.keys())
        print(f"\n{'Model':<12} {'Mean r':>8} {'Median r':>10} {'Pos/Total':>10} "
              f"{'Max r':>8} {'CV R²':>10}")
        print("-" * 62)
        for mname in model_names:
            vals = list(pearson_results[mname].values())
            pos = sum(1 for v in vals if v > 0)
            cv = cv_scores.get(mname, np.nan)
            if isinstance(cv, list):
                cv = np.mean(cv)
            print(f"{mname:<12} {np.mean(vals):>8.4f} {np.median(vals):>10.4f} "
                  f"{pos}/{len(vals):>7} {max(vals):>8.4f} {cv:>10.4f}")

        # Best model per trait — top 15
        print(f"\nTop 15 traits (by best model):")
        best_list = []
        for trait in self.target_traits:
            best_m = max(model_names,
                         key=lambda m: pearson_results[m].get(trait, -999))
            best_r = pearson_results[best_m].get(trait, 0)
            best_list.append((trait, best_m, best_r))
        best_list.sort(key=lambda x: x[2], reverse=True)

        for trait, mname, r in best_list[:15]:
            print(f"  {trait:25s}  {mname:<8s}  r = {r:.4f}")

        n_pos = sum(1 for _, _, r in best_list if r > 0)
        print(f"\n  Positive (best model): {n_pos}/{len(self.target_traits)}")
        print("=" * 70)

    # ──────────────────────────────────────────────────────────────────
    # 12. Main entry point
    # ──────────────────────────────────────────────────────────────────
    def run(self):
        # Load
        pheno_df, traits = self.load_phenotype_data()
        geno_df = self.load_genotype_data(pheno_df[self.accession_col].tolist())

        # Merge & auto-configure
        X, y = self.merge_data()

        # Minimum sample check
        if len(X) < 20:
            logging.error(f"Only {len(X)} samples after merge — too few for ML. "
                          "Check your accession ID matching.")
            sys.exit(1)

        # Train
        results = self.train_ensemble(X, y)

        y_test = results["splits"]["y_test"]
        X_test = results["splits"]["X_test"]
        predictions = results["predictions"]

        # Evaluate
        pearson_results = self.evaluate(y_test, predictions)

        # Uncertainty
        rf_means, rf_stds = self.uncertainty_from_rf(results["models"]["RF"], X_test)
        uncertainty = {"RF": (rf_means, rf_stds)}

        # Feature importance
        imp_dir = os.path.join(self.output_dir, "feature_importance")
        os.makedirs(imp_dir, exist_ok=True)
        rf_imp = self.get_feature_importance(results["models"]["RF"])
        for trait, df in rf_imp.items():
            df.to_csv(os.path.join(imp_dir, f"rf_{trait}.csv"), index=False)

        # Save
        self.save_results(y_test, predictions, pearson_results,
                          results["cv_scores"], uncertainty)

        # Plots
        plot_dir = os.path.join(self.output_dir, "plots")
        # Bar plot for each model
        for mname, pr in pearson_results.items():
            self.plot_pearson_bars(pr, mname, plot_dir)
        # Comparison plot
        self.plot_model_comparison(pearson_results, plot_dir)
        # Scatter plots
        self.plot_predicted_vs_observed(
            y_test, predictions, pearson_results,
            os.path.join(plot_dir, "scatter"),
        )

        logging.info("Pipeline complete")
        return pearson_results


# ═══════════════════════════════════════════════════════════════════════════
# CLI entry point
# ═══════════════════════════════════════════════════════════════════════════
def parse_args():
    """Parse command-line arguments."""
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Adaptive Ensemble ML Pipeline for Genotype-to-Phenotype Prediction.\n"
            "Automatically scales model complexity, PCA components, CV strategy,\n"
            "and ensemble composition based on dataset size."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python adaptive_ensemble_pipeline.py \\\n"
            "      --phenotype data/phenotypes.csv \\\n"
            "      --vcf data/genotypes.vcf\n\n"
            "  python adaptive_ensemble_pipeline.py \\\n"
            "      --phenotype data/rice_traits.csv \\\n"
            "      --vcf data/rice_44k.vcf \\\n"
            "      --accession-col HybID \\\n"
            "      --exclude-cols replicate_id batch \\\n"
            "      --output results_rice/\n"
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

    # ── Optional arguments ────────────────────────────────────────────
    parser.add_argument(
        "--output", default="results",
        help="Output directory for all results (default: results/)",
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
        "--n-components", type=int, default=None,
        help="Number of PCA components (default: auto-determined by regime)",
    )
    parser.add_argument(
        "--test-size", type=float, default=None,
        help="Test set fraction, e.g. 0.2 (default: auto-determined by regime)",
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

    predictor = GenomePhenotypePredictor(
        phenotype_file=args.phenotype,
        vcf_file=args.vcf,
        accession_col=args.accession_col,
        exclude_cols=args.exclude_cols,
        n_components=args.n_components,
        test_size=args.test_size,
        output_dir=args.output,
        remove_outliers=not args.no_outlier_removal,
    )
    predictor.run()
