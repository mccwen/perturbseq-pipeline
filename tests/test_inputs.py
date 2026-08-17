from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import pytest

from perturbseq_pipeline.inputs import prepare_inputs


def _prepare(config: dict, root: Path):
    prepare_inputs(
        config,
        root / "analysis.h5ad",
        root / "expression.mtx",
        root / "genes.tsv",
        root / "barcodes.tsv",
        root / "metadata_out.csv",
    )
    ad = pytest.importorskip("anndata")
    return ad.read_h5ad(root / "analysis.h5ad")


def _write_10x_h5(path: Path, counts: np.ndarray, barcodes: list[str]) -> None:
    h5py = pytest.importorskip("h5py")
    sparse = pytest.importorskip("scipy.sparse")
    feature_by_cell = sparse.csc_matrix(counts)
    names = ["G1", "G2", "guide_1"]
    feature_types = ["Gene Expression", "Gene Expression", "CRISPR Guide Capture"]
    with h5py.File(path, "w") as handle:
        matrix = handle.create_group("matrix")
        matrix.create_dataset("barcodes", data=np.asarray(barcodes, dtype="S"))
        matrix.create_dataset("data", data=feature_by_cell.data.astype(np.int32))
        matrix.create_dataset("indices", data=feature_by_cell.indices.astype(np.int64))
        matrix.create_dataset("indptr", data=feature_by_cell.indptr.astype(np.int64))
        matrix.create_dataset("shape", data=np.asarray(feature_by_cell.shape, dtype=np.int64))
        features = matrix.create_group("features")
        features.create_dataset("id", data=np.asarray(names, dtype="S"))
        features.create_dataset("name", data=np.asarray(names, dtype="S"))
        features.create_dataset("feature_type", data=np.asarray(feature_types, dtype="S"))


def test_prepare_csv_inputs_aligns_guide_cells_and_metadata_defaults():
    with TemporaryDirectory() as directory:
        root = Path(directory)
        pd.DataFrame({"c1": [1], "c2": [2]}, index=["GENE1"]).to_csv(root / "counts.csv")
        pd.DataFrame({"c2": [3], "c1": [4]}, index=["guide_1"]).to_csv(root / "guides.csv")
        pd.DataFrame(
            {
                "cell_barcode": ["c1", "c2"],
                "sample_id": ["S1", "S1"],
                "donor": ["D1", "D1"],
            }
        ).to_csv(root / "metadata.csv", index=False)
        pd.DataFrame(
            {
                "sample_id": ["S1"],
                "donor": ["D1"],
                "batch": ["B1"],
                "condition": ["ctrl"],
                "n_cells": [2],
            }
        ).to_csv(root / "sample_design.csv", index=False)
        config = {
            "input": {
                "mode": "csv",
                "expression_counts": str(root / "counts.csv"),
                "guide_counts": str(root / "guides.csv"),
                "metadata": str(root / "metadata.csv"),
                "sample_design": str(root / "sample_design.csv"),
                "sample_metadata": {"batch": "B1", "condition": "ctrl", "cell_type": "epi"},
            }
        }

        adata = _prepare(config, root)

        sparse = pytest.importorskip("scipy.sparse")
        metadata = pd.read_csv(root / "metadata_out.csv")
        assert sparse.issparse(adata.X)
        assert sparse.issparse(adata.obsm["guide_counts"])
        assert adata.obs_names.tolist() == ["c1", "c2"]
        assert metadata["batch"].tolist() == ["B1", "B1"]


def test_prepare_h5ad_inputs_splits_feature_types_and_uses_obs_metadata():
    ad = pytest.importorskip("anndata")
    with TemporaryDirectory() as directory:
        root = Path(directory)
        adata = ad.AnnData(
            X=np.array([[5, 0, 4], [2, 3, 0]], dtype=np.int64),
            obs=pd.DataFrame(
                {
                    "donor": ["D1", "D2"],
                    "sample_id": ["S1", "S2"],
                    "batch": ["B1", "B1"],
                    "condition": ["ctrl", "ctrl"],
                    "cell_type": ["epi", "epi"],
                },
                index=["c1", "c2"],
            ),
            var=pd.DataFrame(
                {"feature_types": ["Gene Expression", "Gene Expression", "CRISPR Guide Capture"]},
                index=["G1", "G2", "guide_1"],
            ),
        )
        path = root / "input.h5ad"
        adata.write_h5ad(path)
        config = {"input": {"mode": "h5ad", "h5ad": str(path), "metadata": None}}

        result = _prepare(config, root)

        assert result.var_names.tolist() == ["G1", "G2"]
        assert list(result.uns["guide_ids"]) == ["guide_1"]
        assert result.obs["donor"].tolist() == ["D1", "D2"]


def test_prepare_cellranger_h5_reads_feature_barcode_matrix():
    pytest.importorskip("scanpy")
    with TemporaryDirectory() as directory:
        root = Path(directory)
        h5_path = root / "filtered_feature_bc_matrix.h5"
        _write_10x_h5(h5_path, np.array([[5, 2], [0, 3], [4, 0]]), ["AAAC-1", "AAAG-1"])
        metadata = pd.DataFrame(
            {
                "cell_barcode": ["AAAC-1", "AAAG-1"],
                "sample_id": ["S1", "S1"],
                "donor": ["D1", "D1"],
                "batch": ["B1", "B1"],
                "condition": ["vehicle", "vehicle"],
                "cell_type": ["epithelial", "epithelial"],
            }
        )
        metadata.to_csv(root / "metadata.csv", index=False)
        config = {
            "input": {
                "mode": "cellranger_h5",
                "cellranger_h5": str(h5_path),
                "metadata": str(root / "metadata.csv"),
            }
        }

        result = _prepare(config, root)

        assert result[:, "G1"].X.toarray().ravel().tolist() == [5, 2]
        assert result.obsm["guide_counts"].toarray().ravel().tolist() == [4, 0]


def test_manifest_namespaces_repeated_barcodes_across_replicates():
    pytest.importorskip("scanpy")
    with TemporaryDirectory() as directory:
        root = Path(directory)
        for sample_id, guide_count in [("S1", 4), ("S2", 7)]:
            _write_10x_h5(
                root / f"{sample_id}.h5",
                np.array([[5], [2], [guide_count]]),
                ["AAAC-1"],
            )
        pd.DataFrame(
            [
                {
                    "sample_id": "S1",
                    "input_path": "S1.h5",
                    "input_format": "cellranger_h5",
                    "donor": "D1",
                    "batch": "B1",
                    "condition": "vehicle",
                    "cell_type": "epithelial",
                },
                {
                    "sample_id": "S2",
                    "input_path": "S2.h5",
                    "input_format": "cellranger_h5",
                    "donor": "D2",
                    "batch": "B2",
                    "condition": "vehicle",
                    "cell_type": "epithelial",
                },
            ]
        ).to_csv(root / "manifest.csv", index=False)
        config = {"input": {"mode": "manifest", "manifest": str(root / "manifest.csv")}}

        result = _prepare(config, root)

        assert result.obs_names.tolist() == ["S1#AAAC-1", "S2#AAAC-1"]
        assert result.obs["sample_id"].tolist() == ["S1", "S2"]


def test_prepare_inputs_rejects_unresolved_replicate_metadata():
    with TemporaryDirectory() as directory:
        root = Path(directory)
        pd.DataFrame({"c1": [1]}, index=["G1"]).to_csv(root / "counts.csv")
        pd.DataFrame({"c1": [4]}, index=["guide_1"]).to_csv(root / "guides.csv")
        pd.DataFrame(
            {
                "cell_barcode": ["c1"],
                "sample_id": ["S1"],
                "donor": ["unknown"],
                "batch": ["B1"],
                "condition": ["vehicle"],
                "cell_type": ["epithelial"],
            }
        ).to_csv(root / "metadata.csv", index=False)
        config = {
            "input": {
                "mode": "csv",
                "expression_counts": str(root / "counts.csv"),
                "guide_counts": str(root / "guides.csv"),
                "metadata": str(root / "metadata.csv"),
            }
        }

        with pytest.raises(ValueError, match="Required metadata must be resolved"):
            _prepare(config, root)
