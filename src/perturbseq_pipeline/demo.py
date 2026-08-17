from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def make_demo(
    outdir: str | Path,
    n_cells: int = 1200,
    expected_moi: float = 0.30,
    perturbation_response_rate: float = 0.75,
    seed: int = 7,
) -> None:
    """Generate a replicated low-MOI screen with escaped perturbations."""
    rng = np.random.default_rng(seed)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    marker_genes = [
        "TP53",
        "KRAS",
        "MYC",
        "EPCAM",
        "KRT8",
        "MKI67",
        "PTPRC",
        "MT-CO1",
        "RPLP0",
        "ACTB",
        "CDKN1A",
        "BAX",
        "MDM2",
        "DUSP6",
        "ETV4",
        "ETV5",
    ]
    genes = marker_genes + [f"GENE{i:04d}" for i in range(1, 291)]
    sample_rows = []
    for donor in ["D01", "D02", "D03"]:
        for batch in ["B1", "B2"]:
            for condition in ["vehicle", "stimulated"]:
                sample_rows.append(
                    {
                        "sample_id": f"S{len(sample_rows) + 1:02d}",
                        "donor": donor,
                        "batch": batch,
                        "condition": condition,
                    }
                )
    sample_metadata = pd.DataFrame(sample_rows)
    sample_indices = np.arange(n_cells) % len(sample_metadata)
    rng.shuffle(sample_indices)
    sample_ids = sample_metadata.loc[sample_indices, "sample_id"].to_numpy()
    within_sample = pd.Series(sample_ids).groupby(sample_ids).cumcount().to_numpy()
    cells = [
        f"{sample_id}#{index:04d}"
        for sample_id, index in zip(sample_ids, within_sample, strict=True)
    ]
    guides = ["NT_001", "NT_002", "TP53_g1", "KRAS_g1", "MYC_g1"]

    expression = pd.DataFrame(
        rng.negative_binomial(5, 0.35, size=(len(genes), n_cells)),
        index=genes,
        columns=cells,
    )
    guide_counts = pd.DataFrame(0, index=guides, columns=cells)
    perturbation_signatures = {
        "TP53": ["TP53", "CDKN1A", "BAX", "MDM2"] + [f"GENE{i:04d}" for i in range(1, 17)],
        "KRAS": ["KRAS", "DUSP6", "ETV4", "ETV5"] + [f"GENE{i:04d}" for i in range(21, 37)],
        "MYC": ["MYC", "RPLP0"] + [f"GENE{i:04d}" for i in range(41, 59)],
    }
    multiplicities = rng.poisson(expected_moi, size=n_cells)
    for cell, multiplicity in zip(cells, multiplicities):
        if multiplicity == 0:
            continue
        selected = rng.choice(guides, size=min(multiplicity, len(guides)), replace=False)
        for guide in selected:
            guide_counts.loc[guide, cell] = int(rng.integers(3, 25))

        if multiplicity == 1:
            target = str(selected[0]).split("_g", maxsplit=1)[0]
            affected = perturbation_signatures.get(target, [])
            if affected and rng.random() < perturbation_response_rate:
                original = expression.loc[affected, cell].to_numpy(dtype=int)
                expression.loc[affected, cell] = rng.binomial(original, 0.25)

    metadata = sample_metadata.loc[sample_indices].reset_index(drop=True)
    metadata.insert(0, "cell_barcode", cells)
    metadata["cell_type"] = rng.choice(["epithelial", "cycling"], size=n_cells)
    sample_counts = metadata["sample_id"].value_counts()
    sample_metadata["n_cells"] = sample_metadata["sample_id"].map(sample_counts).astype(int)

    expression.to_csv(out / "expression_counts.csv")
    guide_counts.to_csv(out / "guide_counts.csv")
    metadata.to_csv(out / "cell_metadata.csv", index=False)
    sample_metadata.to_csv(out / "sample_metadata.csv", index=False)

    try:
        import anndata as ad
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "AnnData is required to generate the H5AD validation fixture."
        ) from exc

    feature_names = genes + guides
    feature_types = ["Gene Expression"] * len(genes) + ["CRISPR Guide Capture"] * len(guides)
    combined_counts = sparse.hstack(
        [sparse.csr_matrix(expression.T), sparse.csr_matrix(guide_counts.T)],
        format="csr",
    )
    h5ad = ad.AnnData(
        X=combined_counts,
        obs=metadata.set_index("cell_barcode"),
        var=pd.DataFrame({"feature_types": feature_types}, index=feature_names),
    )
    h5ad.obs.index.name = "cell_barcode"
    h5ad.write_h5ad(out / "perturbseq_feature_matrix.h5ad", compression="gzip")
