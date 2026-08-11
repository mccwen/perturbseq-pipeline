import pandas as pd

from perturbseq_pipeline.pseudobulk import aggregate_pseudobulk


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

    pseudobulk, design = aggregate_pseudobulk(
        counts,
        metadata,
        groupby=["donor", "condition", "cell_type", "target_gene"],
        min_cells_per_group=2,
    )

    assert pseudobulk.loc["GENE1"].item() == 4
    assert design.loc[0, "n_cells"] == 2

