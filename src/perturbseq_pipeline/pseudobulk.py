import numpy as np
import pandas as pd
from scipy import sparse

from perturbseq_pipeline.inputs import as_bool


def aggregate_pseudobulk_sparse(
    adata,
    metadata: pd.DataFrame,
    groupby: list[str],
    min_cells_per_group: int = 20,
    eligibility_column: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate sparse cell-by-gene counts before creating a small dense result."""
    if eligibility_column is not None:
        if eligibility_column not in metadata:
            raise ValueError(f"metadata must contain eligibility column: {eligibility_column}")
        metadata = metadata.loc[as_bool(metadata[eligibility_column])].copy()

    meta = metadata.set_index("cell_barcode")
    available = adata.obs_names.intersection(meta.index, sort=False)
    meta = meta.loc[available]
    matrix = sparse.csr_matrix(adata[available, :].X)
    group_keys = meta[groupby].astype(str).agg("|".join, axis=1)
    group_sizes = group_keys.value_counts()
    keep_groups = group_sizes[group_sizes >= min_cells_per_group].index
    keep_cells = group_keys.isin(keep_groups).to_numpy()
    kept_keys = group_keys.loc[keep_cells]
    codes, groups = pd.factorize(kept_keys, sort=True)
    aggregator = sparse.csr_matrix(
        (np.ones(len(codes), dtype=np.int64), (np.arange(len(codes)), codes)),
        shape=(len(codes), len(groups)),
    )
    aggregated = matrix[keep_cells].T @ aggregator
    pseudobulk = pd.DataFrame(
        aggregated.toarray(),
        index=adata.var_names.astype(str),
        columns=groups,
    )
    design = meta.loc[kept_keys.index, groupby].assign(pseudobulk_id=kept_keys.values)
    design = design.drop_duplicates("pseudobulk_id").set_index("pseudobulk_id").loc[groups]
    design.index.name = "pseudobulk_id"
    design = design.reset_index()
    design["n_cells"] = design["pseudobulk_id"].map(group_sizes).astype(int).to_numpy()
    return pseudobulk, design
