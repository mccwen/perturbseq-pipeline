import pandas as pd

from perturbseq_pipeline.qc import calculate_qc_metrics, filter_cells


def test_calculate_qc_metrics_counts_cells_and_mito_fraction():
    counts = pd.DataFrame(
        {
            "cell_a": [10, 0, 5],
            "cell_b": [0, 4, 4],
        },
        index=["GENE1", "MT-CO1", "GENE2"],
    )

    qc = calculate_qc_metrics(counts)

    assert qc.loc[qc["cell_barcode"] == "cell_a", "total_counts"].item() == 15
    assert qc.loc[qc["cell_barcode"] == "cell_b", "detected_genes"].item() == 2
    assert qc.loc[qc["cell_barcode"] == "cell_b", "mito_fraction"].item() == 0.5


def test_filter_cells_adds_pass_qc_flag():
    qc = pd.DataFrame(
        {
            "cell_barcode": ["cell_a", "cell_b"],
            "total_counts": [1000, 100],
            "detected_genes": [500, 50],
            "mito_fraction": [0.05, 0.30],
        }
    )

    filtered = filter_cells(qc, min_counts=500, min_genes=200, max_mito_fraction=0.2)

    assert filtered["pass_qc"].tolist() == [True, False]

