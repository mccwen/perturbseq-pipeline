import numpy as np
import pandas as pd


def expression_pca_sparse(
    adata,
    metadata: pd.DataFrame,
    n_components: int = 10,
    highly_variable_genes: int = 2000,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Return expression PC scores for cells represented in the metadata."""
    try:
        import scanpy as sc
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError("Scanpy is required for sparse expression PCA.") from exc

    cells = adata.obs_names.intersection(metadata["cell_barcode"].astype(str), sort=False)
    if len(cells) < 3:
        raise ValueError("PCA batch diagnosis requires at least three cells")

    analysis = adata[cells, :].copy()
    sc.pp.normalize_total(analysis, target_sum=1e4)
    sc.pp.log1p(analysis)

    n_top = min(highly_variable_genes, analysis.n_vars)
    if n_top < analysis.n_vars:
        sc.pp.highly_variable_genes(
            analysis,
            n_top_genes=n_top,
            flavor="seurat",
            subset=True,
        )

    components = min(n_components, analysis.n_obs - 1, analysis.n_vars - 1)
    if components < 2:
        raise ValueError("PCA batch diagnosis requires at least three cells and genes")
    sc.pp.pca(analysis, n_comps=components, zero_center=True, svd_solver="arpack")

    scores = pd.DataFrame(
        analysis.obsm["X_pca"],
        columns=[f"PC{index + 1}" for index in range(components)],
    )
    scores.insert(0, "cell_barcode", analysis.obs_names.astype(str))
    return scores, np.asarray(analysis.uns["pca"]["variance_ratio"])


def _categorical_r_squared(values: np.ndarray, groups: pd.Series) -> float:
    """Return the fraction of numeric variation associated with a category."""
    grand_mean = float(values.mean())
    total_ss = float(np.square(values - grand_mean).sum())
    if total_ss == 0:
        return 0.0

    between_ss = 0.0
    for indices in groups.groupby(groups).groups.values():
        group_values = values[np.asarray(list(indices), dtype=int)]
        between_ss += len(group_values) * float((group_values.mean() - grand_mean) ** 2)
    return between_ss / total_ss


def contingency_diagnostics(
    metadata: pd.DataFrame,
    covariates: list[str],
    pca_scores: pd.DataFrame | None = None,
    variance_ratio: np.ndarray | None = None,
) -> pd.DataFrame:
    """Summarize target balance and expression association for each covariate."""
    if "target_gene" not in metadata:
        raise ValueError("metadata must contain target_gene")

    rows: list[dict[str, object]] = []
    for covariate in covariates:
        if covariate not in metadata:
            continue

        table = pd.crosstab(metadata[covariate], metadata["target_gene"])
        sparsity = float((table == 0).sum().sum() / table.size) if table.size else 0.0
        weighted_pc_r2 = np.nan
        max_pc_r2 = np.nan

        if pca_scores is not None and variance_ratio is not None:
            aligned = pca_scores.merge(
                metadata[["cell_barcode", covariate]],
                on="cell_barcode",
                validate="one_to_one",
            )
            pc_columns = [column for column in aligned if column.startswith("PC")]
            r_squared = np.array(
                [
                    _categorical_r_squared(aligned[column].to_numpy(), aligned[covariate])
                    for column in pc_columns
                ]
            )
            weights = variance_ratio[: len(r_squared)]
            weighted_pc_r2 = float(np.average(r_squared, weights=weights))
            max_pc_r2 = float(r_squared.max())

        rows.append(
            {
                "covariate": covariate,
                "n_levels": table.shape[0],
                "n_targets": table.shape[1],
                "min_cells_per_level": int(table.sum(axis=1).min()) if len(table) else 0,
                "zero_fraction": round(sparsity, 4),
                "weighted_pc_r2": round(weighted_pc_r2, 4),
                "max_pc_r2": round(max_pc_r2, 4),
                "recommendation": (
                    "include_in_design"
                    if sparsity < 0.5 and (np.isnan(weighted_pc_r2) or weighted_pc_r2 < 0.1)
                    else "inspect_confounding"
                ),
            }
        )
    return pd.DataFrame(rows)


def write_recommendations(diagnostics: pd.DataFrame) -> str:
    """Render expression and balance diagnostics as a short Markdown summary."""
    lines = [
        "# Batch Diagnosis Recommendations",
        "",
        "Use batch correction for visualization only unless the design is demonstrably balanced.",
        "For differential expression, include donor, batch, and condition covariates directly in the model.",
        "The downstream design gate separately blocks exact aliases with donor, guide, or escape status.",
        "",
    ]
    for row in diagnostics.to_dict("records"):
        lines.append(
            f"- `{row['covariate']}`: {row['recommendation']} "
            f"(levels={row['n_levels']}, zero_fraction={row['zero_fraction']}, "
            f"weighted_PC_R2={row['weighted_pc_r2']})"
        )
    return "\n".join(lines) + "\n"
