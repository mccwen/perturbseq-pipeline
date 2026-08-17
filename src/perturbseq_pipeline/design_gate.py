from typing import Any

import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from perturbseq_pipeline.inputs import as_bool

DEFAULT_CONFOUNDERS = [
    "donor",
    "assigned_guide",
    "target_gene",
    "escaped_perturbation",
]


def _connected_components(left: pd.Series, right: pd.Series) -> int:
    """Count disconnected groups in an observed batch-by-covariate design."""
    table = pd.crosstab(left.astype(str), right.astype(str))
    if table.empty:
        return 0
    incidence = sparse.csr_matrix(table.to_numpy())
    graph = sparse.bmat([[None, incidence], [incidence.T, None]], format="csr")
    return int(connected_components(graph, directed=False, return_labels=False))


def evaluate_design_gate(
    metadata: pd.DataFrame,
    batch_column: str = "batch",
    confounders: list[str] | None = None,
    eligibility_column: str | None = "include_primary_analysis",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Block inference when batch-adjusted effects are not identifiable."""
    if batch_column not in metadata:
        raise ValueError(f"metadata must contain batch column: {batch_column}")

    analysis = metadata.copy()
    if eligibility_column:
        if eligibility_column not in analysis:
            raise ValueError(f"metadata must contain eligibility column: {eligibility_column}")
        analysis = analysis.loc[as_bool(analysis[eligibility_column])].copy()
    if analysis.empty:
        raise ValueError("No eligible cells are available for design diagnosis")

    rows: list[dict[str, Any]] = []
    blockers: list[str] = []
    for confounder in confounders or DEFAULT_CONFOUNDERS:
        if confounder not in analysis:
            rows.append(
                {
                    "batch_column": batch_column,
                    "confounder": confounder,
                    "n_batch_levels": analysis[batch_column].nunique(),
                    "n_confounder_levels": 0,
                    "connected_components": 0,
                    "completely_confounded": False,
                    "status": "not_observed",
                }
            )
            continue

        observed = analysis[[batch_column, confounder]].dropna()
        n_batch = observed[batch_column].nunique()
        n_confounder = observed[confounder].nunique()
        components = _connected_components(observed[batch_column], observed[confounder])
        confounded = n_batch > 1 and n_confounder > 1 and components > 1
        if confounded:
            blockers.append(confounder)
        rows.append(
            {
                "batch_column": batch_column,
                "confounder": confounder,
                "n_batch_levels": n_batch,
                "n_confounder_levels": n_confounder,
                "connected_components": components,
                "completely_confounded": confounded,
                "status": "blocked" if confounded else "estimable",
            }
        )

    decision = {
        "status": "blocked" if blockers else "passed",
        "downstream_analysis_allowed": not blockers,
        "batch_column": batch_column,
        "blocking_confounders": blockers,
        "n_eligible_cells": len(analysis),
        "reason": (
            "Batch is exactly aliased with: " + ", ".join(blockers)
            if blockers
            else "No complete batch alias was detected."
        ),
    }
    return pd.DataFrame(rows), decision


def require_design_gate(decision: dict[str, Any]) -> None:
    """Stop aggregation unless the design decision is a pass."""
    if not decision.get("downstream_analysis_allowed", False):
        blockers = ", ".join(decision.get("blocking_confounders", [])) or "unknown"
        raise RuntimeError(
            "Downstream analysis blocked: batch is completely confounded with " + blockers
        )
