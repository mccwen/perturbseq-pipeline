import pandas as pd

from perturbseq_pipeline.guide_calling import assign_guides, summarize_guide_qc


def test_assign_guides_labels_confident_single_and_nt():
    guide_counts = pd.DataFrame(
        {
            "cell_a": [9, 0],
            "cell_b": [1, 7],
        },
        index=["NT_001", "TP53_g1"],
    )
    guide_reference = pd.DataFrame(
        {
            "guide_id": ["NT_001", "TP53_g1"],
            "target_gene": ["NT", "TP53"],
            "is_non_targeting": [True, False],
        }
    )

    calls = assign_guides(guide_counts, guide_reference, min_umi=3, min_fraction=0.6)

    assert calls.loc[calls["cell_barcode"] == "cell_a", "is_non_targeting"].item() is True
    assert calls.loc[calls["cell_barcode"] == "cell_b", "target_gene"].item() == "TP53"
    assert set(calls["guide_assignment"]) == {"single"}


def test_summarize_guide_qc_returns_one_row_per_guide_assignment():
    assignments = pd.DataFrame(
        {
            "cell_barcode": ["c1", "c2"],
            "target_gene": ["TP53", "TP53"],
            "assigned_guide": ["TP53_g1", "TP53_g1"],
            "guide_assignment": ["single", "single"],
            "guide_umi": [10, 6],
            "guide_fraction": [0.8, 0.7],
        }
    )

    summary = summarize_guide_qc(assignments)

    assert summary.loc[0, "n_cells"] == 2
    assert summary.loc[0, "mean_guide_umi"] == 8

