from __future__ import annotations

import pandas as pd


def contingency_diagnostics(metadata: pd.DataFrame, covariates: list[str]) -> pd.DataFrame:
    """Summarize imbalance across covariates and perturbation labels."""
    if "target_gene" not in metadata:
        raise ValueError("metadata must contain target_gene")

    rows: list[dict[str, object]] = []
    for covariate in covariates:
        if covariate not in metadata:
            continue
        table = pd.crosstab(metadata[covariate], metadata["target_gene"])
        min_per_level = int(table.sum(axis=1).min()) if len(table) else 0
        sparsity = float((table == 0).sum().sum() / table.size) if table.size else 0.0
        rows.append(
            {
                "covariate": covariate,
                "n_levels": table.shape[0],
                "n_targets": table.shape[1],
                "min_cells_per_level": min_per_level,
                "zero_fraction": round(sparsity, 4),
                "recommendation": "include_in_design" if sparsity < 0.5 else "inspect_confounding",
            }
        )
    return pd.DataFrame(rows)


def write_recommendations(diagnostics: pd.DataFrame) -> str:
    lines = [
        "# Batch Diagnosis Recommendations",
        "",
        "Use batch correction for visualization only unless the design is demonstrably balanced.",
        "For differential expression, include donor, batch, and condition covariates directly in the model.",
        "",
    ]
    for row in diagnostics.to_dict("records"):
        lines.append(
            f"- `{row['covariate']}`: {row['recommendation']} "
            f"(levels={row['n_levels']}, zero_fraction={row['zero_fraction']})"
        )
    return "\n".join(lines) + "\n"

