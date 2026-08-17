import pandas as pd
from anndata import AnnData
from scipy import sparse

from perturbseq_pipeline.pseudobulk import aggregate_pseudobulk_sparse


def _adata(counts: pd.DataFrame) -> AnnData:
    return AnnData(
        X=sparse.csr_matrix(counts.T),
        obs=pd.DataFrame(index=counts.columns),
        var=pd.DataFrame(index=counts.index),
    )


def test_aggregate_pseudobulk_sums_counts_by_design_group():
    counts = pd.DataFrame(
        {
            "c1": [1, 2],
            "c2": [3, 4],
            "c3": [10, 20],
        },
        index=["GENE1", "GENE2"],
    )
    metadata = pd.DataFrame(
        {
            "cell_barcode": ["c1", "c2", "c3"],
            "donor": ["D1", "D1", "D2"],
            "condition": ["ctrl", "ctrl", "ko"],
            "cell_type": ["epi", "epi", "epi"],
            "target_gene": ["NT", "NT", "TP53"],
        }
    )

    pseudobulk, design = aggregate_pseudobulk_sparse(
        _adata(counts),
        metadata,
        groupby=["donor", "condition", "cell_type", "target_gene"],
        min_cells_per_group=2,
    )

    assert pseudobulk.loc["GENE1"].item() == 4
    assert design.loc[0, "n_cells"] == 2


def test_aggregate_pseudobulk_excludes_low_moi_ineligible_cells():
    counts = pd.DataFrame(
        {"single_1": [1], "single_2": [2], "multiplet": [100]},
        index=["GENE1"],
    )
    metadata = pd.DataFrame(
        {
            "cell_barcode": ["single_1", "single_2", "multiplet"],
            "donor": ["D1", "D1", "D1"],
            "target_gene": ["TP53", "TP53", "unassigned"],
            "include_primary_analysis": [True, True, False],
        }
    )

    pseudobulk, design = aggregate_pseudobulk_sparse(
        _adata(counts),
        metadata,
        groupby=["donor", "target_gene"],
        min_cells_per_group=2,
        eligibility_column="include_primary_analysis",
    )

    assert pseudobulk.loc["GENE1"].item() == 3
    assert design.loc[0, "n_cells"] == 2
