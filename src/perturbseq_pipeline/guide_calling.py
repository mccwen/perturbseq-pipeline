from __future__ import annotations

import numpy as np
import pandas as pd


def assign_guides(
    guide_counts: pd.DataFrame,
    guide_reference: pd.DataFrame,
    min_umi: int = 3,
    min_fraction: float = 0.60,
    max_guides_per_cell: int = 2,
) -> pd.DataFrame:
    """Assign guides from a guides x cells UMI count matrix."""
    guide_to_target = guide_reference.set_index("guide_id")["target_gene"].to_dict()
    nt_guides = set(guide_reference.loc[guide_reference["is_non_targeting"], "guide_id"])
    rows: list[dict[str, object]] = []

    for cell in guide_counts.columns:
        counts = guide_counts[cell].sort_values(ascending=False)
        total = int(counts.sum())
        passing = counts[counts >= min_umi]
        top_guide = str(counts.index[0]) if len(counts) else "unassigned"
        top_umi = int(counts.iloc[0]) if len(counts) else 0
        top_fraction = float(top_umi / total) if total else 0.0
        confident = top_umi >= min_umi and top_fraction >= min_fraction
        n_passing = int((passing > 0).sum())

        if not confident:
            assignment = "unassigned"
        elif n_passing == 1:
            assignment = "single"
        elif n_passing <= max_guides_per_cell:
            assignment = "multi"
        else:
            assignment = "ambiguous"

        rows.append(
            {
                "cell_barcode": cell,
                "assigned_guide": top_guide if confident else "unassigned",
                "target_gene": guide_to_target.get(top_guide, "unassigned") if confident else "unassigned",
                "guide_umi": top_umi,
                "guide_umi_total": total,
                "guide_fraction": top_fraction,
                "n_guides_passing": n_passing,
                "guide_assignment": assignment,
                "is_non_targeting": bool(top_guide in nt_guides and confident),
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

