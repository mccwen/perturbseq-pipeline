from __future__ import annotations

import numpy as np
import pandas as pd


def calculate_qc_metrics(counts: pd.DataFrame, mito_prefix: str = "MT-") -> pd.DataFrame:
    """Return per-cell QC metrics from a genes x cells count matrix."""
    total_counts = counts.sum(axis=0)
    detected_genes = (counts > 0).sum(axis=0)
    mito_genes = counts.index.str.upper().str.startswith(mito_prefix)
    mito_counts = counts.loc[mito_genes].sum(axis=0) if mito_genes.any() else 0
    mito_fraction = np.divide(mito_counts, total_counts, out=np.zeros_like(total_counts, dtype=float), where=total_counts > 0)
    return pd.DataFrame(
        {
            "cell_barcode": counts.columns,
            "total_counts": total_counts.to_numpy(),
            "detected_genes": detected_genes.to_numpy(),
            "mito_fraction": np.asarray(mito_fraction, dtype=float),
        }
    )


def filter_cells(
    qc: pd.DataFrame,
    min_counts: int,
    min_genes: int,
    max_mito_fraction: float,
) -> pd.DataFrame:
    keep = (
        (qc["total_counts"] >= min_counts)
        & (qc["detected_genes"] >= min_genes)
        & (qc["mito_fraction"] <= max_mito_fraction)
    )
    return qc.assign(pass_qc=keep)


def flag_doublets(qc: pd.DataFrame, threshold: float = 0.25) -> pd.DataFrame:
    """Simple placeholder doublet flag based on robust library-size z-score.

    Production workflows can swap this for Scrublet while preserving the output schema.
    """
    total = qc["total_counts"].astype(float)
    median = total.median()
    mad = (total - median).abs().median()
    robust_z = (total - median) / (1.4826 * mad) if mad > 0 else total * 0
    score = 1 / (1 + np.exp(-(robust_z - 3)))
    return qc.assign(doublet_score=score, predicted_doublet=score >= threshold)

