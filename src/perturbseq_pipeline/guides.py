import math

import numpy as np
import pandas as pd
from scipy import sparse


def assign_guides_sparse(
    guide_counts,
    cell_ids: list[str],
    guide_ids: list[str],
    guide_reference: pd.DataFrame,
    min_umi: int = 3,
    min_fraction: float = 0.60,
) -> pd.DataFrame:
    """Assign guides by iterating only over observed entries in a sparse matrix."""
    matrix = sparse.csr_matrix(guide_counts)
    if matrix.shape != (len(cell_ids), len(guide_ids)):
        raise ValueError("Sparse guide matrix dimensions do not match cell and guide identifiers")
    guide_to_target = guide_reference.set_index("guide_id")["target_gene"].to_dict()
    nt_guides = set(guide_reference.loc[guide_reference["is_non_targeting"], "guide_id"])
    rows: list[dict[str, object]] = []
    for row_index, cell in enumerate(cell_ids):
        row = matrix.getrow(row_index)
        order = np.argsort(row.data)[::-1]
        observed_indices = row.indices[order]
        observed_counts = row.data[order]
        total = int(observed_counts.sum())
        top_index = int(observed_indices[0]) if len(observed_indices) else -1
        top_guide = guide_ids[top_index] if top_index >= 0 else "unassigned"
        top_umi = int(observed_counts[0]) if len(observed_counts) else 0
        top_fraction = float(top_umi / total) if total else 0.0
        passing_indices = observed_indices[observed_counts >= min_umi]
        passing_guides = [guide_ids[int(index)] for index in passing_indices]
        passing_targets = sorted(
            {guide_to_target.get(guide, "unknown") for guide in passing_guides}
        )
        n_passing = len(passing_guides)
        confident = top_umi >= min_umi and top_fraction >= min_fraction
        if n_passing > 1:
            assignment = "multiplet"
        elif n_passing == 1 and confident:
            assignment = "single"
        else:
            assignment = "unassigned"
        is_single = assignment == "single"
        rows.append(
            {
                "cell_barcode": cell,
                "assigned_guide": top_guide if is_single else "unassigned",
                "target_gene": (
                    guide_to_target.get(top_guide, "unassigned") if is_single else "unassigned"
                ),
                "guide_umi": top_umi,
                "guide_umi_total": total,
                "guide_fraction": top_fraction,
                "n_guides_passing": n_passing,
                "passing_guides": ";".join(passing_guides),
                "passing_targets": ";".join(passing_targets),
                "guide_assignment": assignment,
                "is_non_targeting": bool(top_guide in nt_guides and is_single),
                "guide_singlet_eligible": is_single,
                "guide_exclusion_reason": "" if is_single else assignment,
            }
        )
    return pd.DataFrame(rows)


def summarize_guide_qc(assignments: pd.DataFrame) -> pd.DataFrame:
    summary = (
        assignments.groupby(["target_gene", "assigned_guide", "guide_assignment"], dropna=False)
        .agg(
            n_cells=("cell_barcode", "count"),
            mean_guide_umi=("guide_umi", "mean"),
            mean_guide_fraction=("guide_fraction", "mean"),
        )
        .reset_index()
    )
    summary["mean_guide_umi"] = np.round(summary["mean_guide_umi"], 3)
    summary["mean_guide_fraction"] = np.round(summary["mean_guide_fraction"], 3)
    return summary.sort_values(["target_gene", "assigned_guide"])


def summarize_low_moi_qc(
    assignments: pd.DataFrame,
    expected_moi: float,
    max_observed_multiplet_fraction: float,
) -> pd.DataFrame:
    """Summarize guide multiplicity under a low-MOI Poisson design.

    ``apparent_moi_from_guide_detection`` is a diagnostic, not an estimate of the
    experimental MOI: guide capture efficiency, selection, and sequencing depth can
    all change the observed guide-positive fraction.
    """
    n_cells = len(assignments)
    n_single = int(assignments["guide_assignment"].eq("single").sum())
    n_multiplet = int(assignments["guide_assignment"].eq("multiplet").sum())
    n_guide_positive = n_single + n_multiplet
    observed_positive_fraction = n_guide_positive / n_cells if n_cells else 0.0
    observed_multiplet_fraction = n_multiplet / n_guide_positive if n_guide_positive else 0.0

    apparent_moi = (
        -math.log1p(-observed_positive_fraction)
        if 0 <= observed_positive_fraction < 1
        else float("nan")
    )
    p_zero = math.exp(-expected_moi)
    p_one = expected_moi * p_zero
    expected_multiplet_fraction = (1 - p_zero - p_one) / (1 - p_zero) if expected_moi > 0 else 0.0

    return pd.DataFrame(
        [
            {
                "perturbation_design": "low_moi_single_guide",
                "expected_moi": expected_moi,
                "n_cells": n_cells,
                "n_guide_negative_or_unassigned": n_cells - n_guide_positive,
                "n_single_guide": n_single,
                "n_guide_multiplet": n_multiplet,
                "observed_guide_positive_fraction": round(observed_positive_fraction, 4),
                "observed_multiplet_fraction_among_guide_positive": round(
                    observed_multiplet_fraction, 4
                ),
                "poisson_expected_multiplet_fraction_at_target_moi": round(
                    expected_multiplet_fraction, 4
                ),
                "apparent_moi_from_guide_detection": round(apparent_moi, 4),
                "max_observed_multiplet_fraction": max_observed_multiplet_fraction,
                "passes_low_moi_qc": observed_multiplet_fraction <= max_observed_multiplet_fraction,
                "apparent_moi_caveat": (
                    "Observed guide detection is affected by capture, selection, and sequencing depth."
                ),
            }
        ]
    )
