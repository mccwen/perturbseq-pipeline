import pandas as pd
from anndata import AnnData
from scipy import sparse

from perturbseq_pipeline.qc import calculate_qc_metrics, call_cells, filter_cells


def test_calculate_qc_metrics_counts_cells_and_mito_fraction():
    counts = pd.DataFrame(
        {
            "cell_a": [10, 0, 5],
            "cell_b": [0, 4, 4],
        },
        index=["GENE1", "MT-CO1", "GENE2"],
    )

    adata = AnnData(
        X=sparse.csr_matrix(counts.T),
        obs=pd.DataFrame(index=counts.columns),
        var=pd.DataFrame(index=counts.index),
    )
    qc = calculate_qc_metrics(adata)

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


def test_fixed_threshold_cell_calling_is_part_of_qc_gate():
    qc = pd.DataFrame(
        {
            "cell_barcode": ["cell", "empty_droplet"],
            "total_counts": [1000, 20],
            "detected_genes": [500, 10],
            "mito_fraction": [0.05, 0.01],
        }
    )

    called = call_cells(qc, method="fixed_threshold", min_counts=100, min_genes=50)
    filtered = filter_cells(called, min_counts=10, min_genes=5, max_mito_fraction=0.2)

    assert called["called_cell"].tolist() == [True, False]
    assert filtered["pass_qc"].tolist() == [True, False]
