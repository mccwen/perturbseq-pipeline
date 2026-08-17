import numpy as np
import pandas as pd
from scipy import sparse


def calculate_qc_metrics(counts, mito_prefix: str = "MT-") -> pd.DataFrame:
    """Return per-cell QC metrics without densifying sparse analysis matrices."""
    matrix = sparse.csr_matrix(counts.X)
    genes = counts.var_names.astype(str)
    cells = counts.obs_names.astype(str)
    total_counts = np.asarray(matrix.sum(axis=1)).ravel()
    detected_genes = matrix.getnnz(axis=1)
    mito_genes = genes.str.upper().str.startswith(mito_prefix)
    mito_counts = (
        np.asarray(matrix[:, mito_genes].sum(axis=1)).ravel()
        if mito_genes.any()
        else np.zeros(matrix.shape[0])
    )
    mito_fraction = np.divide(
        mito_counts,
        total_counts,
        out=np.zeros_like(total_counts, dtype=float),
        where=total_counts > 0,
    )
    return pd.DataFrame(
        {
            "cell_barcode": cells,
            "total_counts": total_counts,
            "detected_genes": detected_genes,
            "mito_fraction": np.asarray(mito_fraction, dtype=float),
        }
    )


def call_cells(
    qc: pd.DataFrame,
    method: str = "filtered_input",
    min_counts: int = 100,
    min_genes: int = 50,
) -> pd.DataFrame:
    """Record or apply the configured cell-calling boundary.

    ``filtered_input`` trusts Cell Ranger's filtered matrix or a prefiltered h5ad.
    ``fixed_threshold`` is provided for already-demultiplexed custom matrices and makes
    the calling rule fully auditable; raw droplet matrices should normally be called by
    Cell Ranger before entering this pipeline.
    """
    if method == "filtered_input":
        called = pd.Series(True, index=qc.index)
    elif method == "fixed_threshold":
        called = (qc["total_counts"] >= min_counts) & (qc["detected_genes"] >= min_genes)
    else:
        raise ValueError("cell calling method must be filtered_input or fixed_threshold")
    return qc.assign(
        called_cell=called.to_numpy(),
        cell_calling_method=method,
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
    if "called_cell" in qc:
        keep &= qc["called_cell"]
    return qc.assign(pass_qc=keep)


def score_doublets_scrublet(
    counts,
    qc: pd.DataFrame,
    expected_doublet_rate: float = 0.06,
    threshold: float | None = None,
) -> pd.DataFrame:
    """Attach Scrublet doublet scores and calls to a QC table."""
    try:
        import scrublet as scr
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "Scrublet is required for doublet detection. Install the Conda environment."
        ) from exc

    matrix = sparse.csr_matrix(counts.X)
    cells = counts.obs_names.astype(str)
    n_cells, n_genes = matrix.shape
    n_prin_comps = min(30, n_cells - 1, n_genes - 1)
    if n_prin_comps < 2:
        raise ValueError("Scrublet requires at least three cells and three genes")

    scrublet = scr.Scrublet(
        matrix,
        expected_doublet_rate=expected_doublet_rate,
        random_state=0,
    )
    scores, automatic_calls = scrublet.scrub_doublets(
        min_counts=2,
        min_cells=3,
        min_gene_variability_pctl=85,
        n_prin_comps=n_prin_comps,
        verbose=False,
    )
    calls = scores >= threshold if threshold is not None else automatic_calls
    score_table = pd.DataFrame(
        {
            "cell_barcode": cells,
            "doublet_score": scores,
            "predicted_doublet": calls,
        }
    )
    return qc.merge(score_table, on="cell_barcode", how="left", validate="one_to_one")
