from __future__ import annotations

import pandas as pd


def aggregate_pseudobulk(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    groupby: list[str],
    min_cells_per_group: int = 20,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate genes x cells counts into genes x pseudobulk samples."""
    meta = metadata.set_index("cell_barcode")
    available_cells = [cell for cell in counts.columns if cell in meta.index]
    counts = counts[available_cells]
    meta = meta.loc[available_cells]

    group_keys = meta[groupby].astype(str).agg("|".join, axis=1)
    group_sizes = group_keys.value_counts()
    keep_groups = group_sizes[group_sizes >= min_cells_per_group].index
    kept_cells = group_keys[group_keys.isin(keep_groups)]

    pseudobulk = counts[kept_cells.index].T.groupby(kept_cells).sum().T
    design = meta.loc[kept_cells.index, groupby].assign(pseudobulk_id=kept_cells.values)
    design = design.drop_duplicates("pseudobulk_id").set_index("pseudobulk_id").loc[pseudobulk.columns]
    design.index.name = "pseudobulk_id"
    design = design.reset_index()
    design["n_cells"] = design["pseudobulk_id"].map(group_sizes).astype(int).to_numpy()
    return pseudobulk, design
