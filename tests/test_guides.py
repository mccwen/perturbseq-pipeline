import pandas as pd
from scipy import sparse

from perturbseq_pipeline.guides import (
    assign_guides_sparse,
    summarize_guide_qc,
    summarize_low_moi_qc,
)


def _assign(guide_counts: pd.DataFrame, guide_reference: pd.DataFrame) -> pd.DataFrame:
    return assign_guides_sparse(
        sparse.csr_matrix(guide_counts.T),
        guide_counts.columns.tolist(),
        guide_counts.index.tolist(),
        guide_reference,
        min_umi=3,
        min_fraction=0.6,
    )


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

    calls = _assign(guide_counts, guide_reference)

    assert calls.loc[calls["cell_barcode"] == "cell_a", "is_non_targeting"].item() is True
    assert calls.loc[calls["cell_barcode"] == "cell_b", "target_gene"].item() == "TP53"
    assert set(calls["guide_assignment"]) == {"single"}
    assert calls["guide_singlet_eligible"].all()


def test_low_moi_guide_multiplet_is_not_assigned_to_one_target():
    guide_counts = pd.DataFrame(
        {"cell_multiplet": [8, 7, 0]},
        index=["NT_001", "TP53_g1", "KRAS_g1"],
    )
    guide_reference = pd.DataFrame(
        {
            "guide_id": ["NT_001", "TP53_g1", "KRAS_g1"],
            "target_gene": ["NT", "TP53", "KRAS"],
            "is_non_targeting": [True, False, False],
        }
    )

    call = _assign(guide_counts, guide_reference).iloc[0]

    assert call["guide_assignment"] == "multiplet"
    assert call["target_gene"] == "unassigned"
    assert call["passing_targets"] == "NT;TP53"
    assert bool(call["guide_singlet_eligible"]) is False


def test_low_moi_qc_reports_collision_rate_and_apparent_moi():
    assignments = pd.DataFrame(
        {"guide_assignment": ["single"] * 20 + ["multiplet"] * 2 + ["unassigned"] * 78}
    )

    summary = summarize_low_moi_qc(
        assignments,
        expected_moi=0.3,
        max_observed_multiplet_fraction=0.2,
    ).iloc[0]

    assert summary["n_single_guide"] == 20
    assert summary["n_guide_multiplet"] == 2
    assert summary["observed_multiplet_fraction_among_guide_positive"] == 0.0909
    assert bool(summary["passes_low_moi_qc"]) is True


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
