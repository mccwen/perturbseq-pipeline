import json

import anndata as ad
import pandas as pd
from scipy import sparse

from perturbseq_pipeline.provenance import file_fingerprint, write_run_manifest


def test_file_fingerprint_records_size_and_sha256(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text("perturb-seq\n")

    fingerprint = file_fingerprint(source, tmp_path)

    assert fingerprint["path"] == "./input.txt"
    assert fingerprint["size_bytes"] == 12
    assert len(fingerprint["sha256"]) == 64


def test_run_manifest_records_inputs_and_sparse_shape(tmp_path):
    root = tmp_path
    (root / "config").mkdir()
    (root / "demo").mkdir()
    for name in ["expression.csv", "guides.csv", "metadata.csv", "reference.csv", "sets.gmt"]:
        (root / "demo" / name).write_text(f"fixture:{name}\n")
    config = root / "config" / "study.yaml"
    config.write_text(
        """input:
  mode: csv
  expression_counts: demo/expression.csv
  guide_counts: demo/guides.csv
  metadata: demo/metadata.csv
guide_assignment:
  reference: demo/reference.csv
pathway_analysis:
  gene_sets: demo/sets.gmt
"""
    )
    adata = ad.AnnData(
        X=sparse.csr_matrix([[1, 0], [0, 2]]),
        obs=pd.DataFrame(index=["c1", "c2"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    adata.layers["counts"] = adata.X.copy()
    adata.obsm["guide_counts"] = sparse.csr_matrix([[3], [0]])
    adata.uns["guide_ids"] = ["NT_1"]
    adata.uns["pipeline_schema_version"] = "1.0"
    analysis = root / "analysis.h5ad"
    adata.write_h5ad(analysis)

    output = root / "results" / "run_metadata.json"
    manifest = write_run_manifest(config, analysis, output)

    saved = json.loads(output.read_text())
    assert saved["schema_version"] == "1.1"
    assert saved["configuration"]["path"] == "./config/study.yaml"
    assert saved["analysis_object"]["path"] == "./analysis.h5ad"
    assert all(item["path"].startswith("./demo/") for item in saved["source_inputs"])
    assert "created_utc" not in saved
    assert "execution" not in saved
    assert manifest["analysis_object"]["n_cells"] == 2
    assert saved["analysis_object"]["expression_nonzero"] == 2
    assert len(saved["source_inputs"]) == 5
