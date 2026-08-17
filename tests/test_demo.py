from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from perturbseq_pipeline.demo import make_demo
from perturbseq_pipeline.inputs import prepare_inputs


def test_demo_uses_sample_level_fully_crossed_replicates():
    with TemporaryDirectory() as directory:
        make_demo(directory, n_cells=120)
        samples = pd.read_csv(f"{directory}/sample_metadata.csv")
        cells = pd.read_csv(f"{directory}/cell_metadata.csv")

        assert len(samples) == 12
        assert len(samples[["donor", "batch", "condition"]].drop_duplicates()) == 12
        observed = cells.merge(
            samples[["sample_id", "donor", "batch", "condition"]],
            on="sample_id",
            suffixes=("", "_sample"),
            validate="many_to_one",
        )
        for column in ["donor", "batch", "condition"]:
            assert observed[column].equals(observed[f"{column}_sample"])
        barcode_samples = cells["cell_barcode"].str.split("#", n=1).str[0]
        assert barcode_samples.equals(cells["sample_id"])


def _normalize(config: dict, output: Path) -> ad.AnnData:
    prepare_inputs(
        config,
        output / "analysis.h5ad",
        output / "expression.mtx",
        output / "genes.tsv",
        output / "barcodes.tsv",
        output / "metadata.csv",
    )
    return ad.read_h5ad(output / "analysis.h5ad")


def test_demo_csv_and_h5ad_normalize_to_identical_sparse_objects():
    with TemporaryDirectory() as directory:
        root = Path(directory)
        make_demo(root, n_cells=120)
        csv_result = _normalize(
            {
                "input": {
                    "mode": "csv",
                    "expression_counts": root / "expression_counts.csv",
                    "guide_counts": root / "guide_counts.csv",
                    "metadata": root / "cell_metadata.csv",
                    "sample_design": root / "sample_metadata.csv",
                }
            },
            root / "csv",
        )
        h5ad_result = _normalize(
            {
                "input": {
                    "mode": "h5ad",
                    "h5ad": root / "perturbseq_feature_matrix.h5ad",
                    "sample_design": root / "sample_metadata.csv",
                }
            },
            root / "h5ad",
        )

        assert sparse.issparse(csv_result.X)
        assert sparse.issparse(h5ad_result.X)
        assert np.array_equal(csv_result.X.toarray(), h5ad_result.X.toarray())
        assert np.array_equal(
            csv_result.obsm["guide_counts"].toarray(),
            h5ad_result.obsm["guide_counts"].toarray(),
        )
        assert csv_result.var_names.equals(h5ad_result.var_names)
        assert csv_result.obs_names.equals(h5ad_result.obs_names)
        pd.testing.assert_frame_equal(csv_result.obs, h5ad_result.obs, check_dtype=False)
