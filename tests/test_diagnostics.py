import pandas as pd
from anndata import AnnData
from scipy import sparse

from perturbseq_pipeline.expression_diagnostics import (
    contingency_diagnostics,
    expression_pca_sparse,
)
from perturbseq_pipeline.guide_diagnostics import guide_batch_diagnostics


def test_expression_pca_and_batch_association_are_reported():
    counts = pd.DataFrame(
        {
            "c1": [10, 0, 1],
            "c2": [9, 0, 2],
            "c3": [0, 10, 1],
            "c4": [0, 9, 2],
        },
        index=["G1", "G2", "G3"],
    )
    metadata = pd.DataFrame(
        {
            "cell_barcode": ["c1", "c2", "c3", "c4"],
            "batch": ["B1", "B1", "B2", "B2"],
            "target_gene": ["NT", "TP53", "NT", "TP53"],
        }
    )

    adata = AnnData(
        X=sparse.csr_matrix(counts.T),
        obs=pd.DataFrame(index=counts.columns),
        var=pd.DataFrame(index=counts.index),
    )
    scores, variance = expression_pca_sparse(
        adata,
        metadata,
        n_components=2,
        highly_variable_genes=3,
    )
    diagnostics = contingency_diagnostics(metadata, ["batch"], scores, variance)

    assert scores.shape == (4, 3)
    assert diagnostics.loc[0, "max_pc_r2"] > 0.9


def test_guide_batch_diagnostics_use_samples_and_flag_large_batch_effect():
    metadata_rows = []
    assignment_rows = []
    for batch in ["B1", "B2"]:
        for replicate in range(5):
            sample_id = f"{batch}_S{replicate}"
            for cell_index in range(10):
                barcode = f"{sample_id}_C{cell_index}"
                is_single = batch == "B1"
                metadata_rows.append(
                    {"cell_barcode": barcode, "sample_id": sample_id, "batch": batch}
                )
                assignment_rows.append(
                    {
                        "cell_barcode": barcode,
                        "assigned_guide": "TP53_g1" if is_single else "unassigned",
                        "guide_assignment": "single" if is_single else "unassigned",
                        "guide_umi_total": 10 if is_single else 0,
                        "guide_fraction": 0.9 if is_single else 0.0,
                    }
                )

    by_sample, tests, composition = guide_batch_diagnostics(
        pd.DataFrame(assignment_rows),
        pd.DataFrame(metadata_rows),
        min_samples_per_batch=2,
    )

    assert len(by_sample) == 10
    positive_test = tests.loc[tests["metric"].eq("guide_positive_rate")].iloc[0]
    assert positive_test["status"] == "tested"
    assert positive_test["min_samples_per_batch"] == 5
    assert positive_test["fdr"] < 0.05
    assert bool(positive_test["batch_effect_detected"]) is True
    assert set(composition["batch"]) == {"B1", "B2"}


def test_guide_batch_diagnostics_report_insufficient_replication():
    metadata = pd.DataFrame(
        {
            "cell_barcode": ["c1", "c2"],
            "sample_id": ["S1", "S2"],
            "batch": ["B1", "B2"],
        }
    )
    assignments = pd.DataFrame(
        {
            "cell_barcode": ["c1", "c2"],
            "assigned_guide": ["NT_g1", "NT_g1"],
            "guide_assignment": ["single", "single"],
            "guide_umi_total": [8, 9],
            "guide_fraction": [0.8, 0.9],
        }
    )

    _, tests, _ = guide_batch_diagnostics(assignments, metadata)

    assert tests["status"].eq("insufficient_replication").all()
    assert tests["fdr"].isna().all()
    assert not tests["batch_effect_detected"].any()
