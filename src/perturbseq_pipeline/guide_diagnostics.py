from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, kruskal

GUIDE_METRICS = [
    "guide_positive_rate",
    "singlet_rate",
    "unassigned_rate",
    "multiplet_rate_among_guide_positive",
    "mean_guide_umi_total",
    "median_guide_umi_positive",
    "median_top_guide_fraction",
]


def _summarize_samples(
    frame: pd.DataFrame,
    sample_column: str,
    batch_column: str,
) -> pd.DataFrame:
    """Calculate guide-capture metrics for each biological sample."""
    rows: list[dict[str, Any]] = []
    for (sample_id, batch), sample in frame.groupby([sample_column, batch_column], sort=True):
        assignment = sample["guide_assignment"]
        positive = assignment.isin(["single", "multiplet"])
        detected = sample["guide_umi_total"].gt(0)
        n_positive = int(positive.sum())
        n_multiplet = int(assignment.eq("multiplet").sum())

        rows.append(
            {
                sample_column: sample_id,
                batch_column: batch,
                "n_cells": len(sample),
                "n_guide_positive": n_positive,
                "n_single": int(assignment.eq("single").sum()),
                "n_multiplet": n_multiplet,
                "n_unassigned": int(assignment.eq("unassigned").sum()),
                "guide_positive_rate": float(positive.mean()),
                "singlet_rate": float(assignment.eq("single").mean()),
                "unassigned_rate": float(assignment.eq("unassigned").mean()),
                "multiplet_rate_among_guide_positive": (
                    float(n_multiplet / n_positive) if n_positive else np.nan
                ),
                "mean_guide_umi_total": float(sample["guide_umi_total"].mean()),
                "median_guide_umi_positive": (
                    float(sample.loc[positive, "guide_umi_total"].median())
                    if n_positive
                    else np.nan
                ),
                "median_top_guide_fraction": (
                    float(sample.loc[detected, "guide_fraction"].median())
                    if detected.any()
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def _batch_test(
    values: pd.DataFrame,
    metric: str,
    scope: str,
    batch_column: str,
    min_samples_per_batch: int,
) -> dict[str, Any]:
    """Run one replicate-level Kruskal-Wallis batch comparison."""
    groups = [
        group[metric].dropna().to_numpy(dtype=float)
        for _, group in values.groupby(batch_column, sort=True)
    ]
    n_batches = len(groups)
    min_replicates = min((len(group) for group in groups), default=0)
    row: dict[str, Any] = {
        "scope": scope,
        "metric": metric,
        "n_batches": n_batches,
        "n_samples": int(sum(len(group) for group in groups)),
        "min_samples_per_batch": min_replicates,
        "statistic": np.nan,
        "p_value": np.nan,
        "effect_size_epsilon_squared": np.nan,
        "status": "insufficient_replication",
    }
    if n_batches < 2 or min_replicates < min_samples_per_batch:
        return row

    combined = np.concatenate(groups)
    statistic, p_value = (0.0, 1.0) if np.allclose(combined, combined[0]) else kruskal(*groups)
    denominator = len(combined) - n_batches
    effect_size = max(0.0, (statistic - n_batches + 1) / denominator)
    row.update(
        statistic=float(statistic),
        p_value=float(p_value),
        effect_size_epsilon_squared=float(effect_size),
        status="tested",
    )
    return row


def guide_batch_diagnostics(
    assignments: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_column: str = "sample_id",
    batch_column: str = "batch",
    min_samples_per_batch: int = 2,
    fdr_threshold: float = 0.05,
    effect_size_threshold: float = 0.10,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Test sample-level guide-capture and assignment metrics across batches."""
    required_assignments = {
        "cell_barcode",
        "assigned_guide",
        "guide_assignment",
        "guide_umi_total",
        "guide_fraction",
    }
    required_metadata = {"cell_barcode", sample_column, batch_column}
    missing_assignments = required_assignments - set(assignments)
    missing_metadata = required_metadata - set(metadata)
    if missing_assignments:
        raise ValueError(f"Guide assignments are missing columns: {sorted(missing_assignments)}")
    if missing_metadata:
        raise ValueError(f"Metadata is missing columns: {sorted(missing_metadata)}")
    if min_samples_per_batch < 2:
        raise ValueError("min_samples_per_batch must be at least 2")
    if metadata.groupby(sample_column)[batch_column].nunique().gt(1).any():
        raise ValueError("Each sample_id must belong to exactly one batch")

    frame = assignments.merge(
        metadata[["cell_barcode", sample_column, batch_column]],
        on="cell_barcode",
        how="inner",
        validate="one_to_one",
    )
    if len(frame) != len(assignments):
        raise ValueError("Every guide assignment must have matching sample and batch metadata")

    by_sample = _summarize_samples(frame, sample_column, batch_column)
    test_rows = [
        _batch_test(by_sample, metric, "capture_metric", batch_column, min_samples_per_batch)
        for metric in GUIDE_METRICS
    ]

    singles = frame.loc[frame["guide_assignment"].eq("single")].copy()
    guides = sorted(singles["assigned_guide"].dropna().unique())
    sample_index = by_sample[[sample_column, batch_column]].drop_duplicates()
    sample_single_totals = singles.groupby(sample_column).size()
    for guide in guides:
        guide_counts = (
            singles.loc[singles["assigned_guide"].eq(guide)].groupby(sample_column).size()
        )
        guide_values = sample_index.copy()
        guide_values["guide_recovery_rate"] = guide_values[sample_column].map(guide_counts).fillna(
            0
        ) / guide_values[sample_column].map(sample_single_totals).replace(0, np.nan)
        test_rows.append(
            _batch_test(
                guide_values,
                "guide_recovery_rate",
                f"guide:{guide}",
                batch_column,
                min_samples_per_batch,
            )
        )

    tests = pd.DataFrame(test_rows)
    tests["fdr"] = np.nan
    tested = tests["p_value"].notna()
    tests.loc[tested, "fdr"] = false_discovery_control(
        tests.loc[tested, "p_value"].to_numpy(),
        method="bh",
    )
    tests["batch_effect_detected"] = (
        tests["status"].eq("tested")
        & tests["fdr"].le(fdr_threshold)
        & tests["effect_size_epsilon_squared"].ge(effect_size_threshold)
    )
    tests["fdr_threshold"] = fdr_threshold
    tests["effect_size_threshold"] = effect_size_threshold

    batch_totals = singles.groupby(batch_column).size().rename("n_single_cells")
    composition = (
        singles.groupby([batch_column, "assigned_guide"])
        .size()
        .rename("n_assigned_cells")
        .reindex(
            pd.MultiIndex.from_product(
                [sorted(frame[batch_column].unique()), guides],
                names=[batch_column, "assigned_guide"],
            ),
            fill_value=0,
        )
        .reset_index()
    )
    composition["n_single_cells"] = (
        composition[batch_column].map(batch_totals).fillna(0).astype(int)
    )
    composition["guide_fraction_among_singlets"] = np.divide(
        composition["n_assigned_cells"],
        composition["n_single_cells"],
        out=np.zeros(len(composition), dtype=float),
        where=composition["n_single_cells"].to_numpy() > 0,
    )
    return by_sample, tests, composition
