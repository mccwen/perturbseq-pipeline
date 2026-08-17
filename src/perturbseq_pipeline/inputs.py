from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from scipy.io import mmwrite

REQUIRED_METADATA = [
    "cell_barcode",
    "sample_id",
    "donor",
    "batch",
    "condition",
    "cell_type",
]
INVALID_METADATA_VALUES = {"", "unknown", "unassigned", "na", "nan", "none", "null"}
MANIFEST_REQUIRED = {"sample_id", "input_path", "input_format", "donor", "batch", "condition"}
SAMPLE_DESIGN_REQUIRED = {"sample_id", "donor", "batch", "condition"}


def read_yaml(path: str | Path) -> dict[str, Any]:
    with Path(path).open() as handle:
        return yaml.safe_load(handle)


def ensure_parent(path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    return out


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin(["true", "1", "yes"])


def resolve_path(path: str | Path, base: Path) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else (base / candidate).resolve()


def configured_input_paths(config: dict[str, Any], root: str | Path = ".") -> list[Path]:
    """Resolve primary matrix and metadata inputs for Snakemake and provenance."""
    root = Path(root).resolve()
    input_config = config["input"]
    mode = input_config["mode"]
    if mode == "csv":
        sources = [
            resolve_path(input_config["expression_counts"], root),
            resolve_path(input_config["guide_counts"], root),
            resolve_path(input_config["metadata"], root),
        ]
    elif mode in {"cellranger_h5", "h5ad"}:
        sources = [resolve_path(input_config[mode], root)]
        if input_config.get("metadata"):
            sources.append(resolve_path(input_config["metadata"], root))
    elif mode == "manifest":
        manifest_path = resolve_path(input_config["manifest"], root)
        sources = [manifest_path]
        manifest = pd.read_csv(manifest_path, keep_default_na=False)
        for row in manifest.to_dict("records"):
            for column in ("input_path", "metadata_path"):
                value = str(row.get(column, "")).strip()
                if value:
                    sources.append(resolve_path(value, manifest_path.parent))
    else:
        raise ValueError("input.mode must be csv, cellranger_h5, h5ad, or manifest")
    if input_config.get("sample_design"):
        sources.append(resolve_path(input_config["sample_design"], root))
    return list(dict.fromkeys(sources))


def read_counts_csv(path: str | Path) -> pd.DataFrame:
    """Read a features-by-cells CSV fixture."""
    counts = pd.read_csv(path, index_col=0)
    counts.index = counts.index.astype(str)
    counts.columns = counts.columns.astype(str)
    return counts


def write_table(df: pd.DataFrame, path: str | Path) -> None:
    df.to_csv(ensure_parent(path), index=False)


def write_matrix(df: pd.DataFrame, path: str | Path) -> None:
    df.to_csv(ensure_parent(path))


def read_guides(path: str | Path) -> pd.DataFrame:
    guides = pd.read_csv(path)
    required = {"guide_id", "target_gene", "is_non_targeting"}
    missing = required - set(guides.columns)
    if missing:
        raise ValueError(f"guide reference is missing columns: {sorted(missing)}")
    guides["is_non_targeting"] = as_bool(guides["is_non_targeting"])
    return guides


def _as_count_csr(matrix: Any, label: str) -> sparse.csr_matrix:
    result = sparse.csr_matrix(matrix)
    if result.shape[0] == 0 or result.shape[1] == 0:
        raise ValueError(f"{label} is empty")
    if not np.issubdtype(result.data.dtype, np.number) or np.isnan(result.data).any():
        raise ValueError(f"{label} must contain numeric counts without missing values")
    if (result.data < 0).any() or not np.allclose(result.data, np.round(result.data)):
        raise ValueError(f"{label} must contain non-negative integer counts")
    result.data = np.round(result.data).astype(np.int64)
    result.eliminate_zeros()
    return result


def _prepare_metadata(
    metadata: pd.DataFrame | None,
    cells: pd.Index,
    defaults: dict[str, Any] | None,
) -> pd.DataFrame:
    defaults = defaults or {}
    supplied_metadata = metadata is not None
    metadata = pd.DataFrame({"cell_barcode": cells}) if metadata is None else metadata.copy()
    if "cell_barcode" not in metadata:
        raise ValueError("Cell metadata must contain cell_barcode")
    if metadata["cell_barcode"].duplicated().any():
        raise ValueError("Cell metadata contains duplicate cell_barcode values")

    metadata["cell_barcode"] = metadata["cell_barcode"].astype(str)
    if supplied_metadata:
        missing_cells = cells.difference(pd.Index(metadata["cell_barcode"]))
        if len(missing_cells) > 0:
            preview = ", ".join(missing_cells[:5])
            raise ValueError(f"Cell metadata is missing {len(missing_cells)} barcodes: {preview}")

    metadata = metadata.set_index("cell_barcode").reindex(cells)
    for column in REQUIRED_METADATA[1:]:
        fallback = defaults.get(column, "unknown")
        fallback = "unknown" if fallback is None else fallback
        if column not in metadata:
            metadata[column] = fallback
        else:
            metadata[column] = metadata[column].astype(object).fillna(fallback)
    invalid_columns = []
    for column in REQUIRED_METADATA[1:]:
        values = metadata[column].astype(str).str.strip()
        invalid = values.str.lower().isin(INVALID_METADATA_VALUES)
        if invalid.any():
            invalid_columns.append(f"{column} ({int(invalid.sum())} cells)")
        metadata[column] = values
    if invalid_columns:
        raise ValueError(
            "Required metadata must be resolved before analysis: " + ", ".join(invalid_columns)
        )
    metadata.index.name = "cell_barcode"
    return metadata


def _validate_sample_design(metadata: pd.DataFrame, path: str | Path) -> None:
    design = pd.read_csv(path)
    missing = SAMPLE_DESIGN_REQUIRED - set(design)
    if missing:
        raise ValueError(f"Sample design is missing columns: {sorted(missing)}")
    if design.empty or design["sample_id"].duplicated().any():
        raise ValueError("Sample design must contain one row per sample_id")

    columns = ["sample_id", "donor", "batch", "condition"]
    observed = metadata.reset_index()[columns]
    inconsistent = observed.groupby("sample_id")[["donor", "batch", "condition"]].nunique()
    if inconsistent.gt(1).any().any():
        raise ValueError("Per-cell metadata assigns a sample to multiple design values")
    observed = observed.drop_duplicates("sample_id").set_index("sample_id").sort_index()
    expected = design[columns].astype(str).set_index("sample_id").sort_index()
    if not observed.astype(str).equals(expected):
        raise ValueError(
            "Sample design does not match per-cell donor, batch, or condition metadata"
        )
    if "n_cells" in design:
        observed_counts = metadata["sample_id"].value_counts().sort_index()
        expected_counts = design.set_index("sample_id")["n_cells"].astype(int).sort_index()
        if not observed_counts.equals(expected_counts):
            raise ValueError("Sample design n_cells does not match per-cell metadata")


def _split_feature_types(adata: Any) -> Any:
    """Create a sparse RNA AnnData with sparse guide counts in ``obsm``."""
    adata.var_names_make_unique()
    feature_column = next(
        (column for column in ("feature_types", "feature_type") if column in adata.var),
        None,
    )
    if feature_column is None:
        raise ValueError(
            "AnnData must contain var['feature_types'] (or var['feature_type']) with both "
            "Gene Expression and CRISPR Guide Capture features"
        )
    feature_types = adata.var[feature_column].astype(str)
    gene_mask = feature_types.eq("Gene Expression").to_numpy()
    guide_mask = feature_types.str.contains(
        "CRISPR|Guide Capture", case=False, regex=True
    ).to_numpy()
    if not gene_mask.any():
        raise ValueError("Input contains no Gene Expression features")
    if not guide_mask.any():
        raise ValueError("Input contains no CRISPR Guide Capture features")

    result = ad.AnnData(
        X=_as_count_csr(adata[:, gene_mask].X, "expression count matrix"),
        obs=adata.obs.copy(),
        var=adata.var.loc[gene_mask].copy(),
    )
    result.obsm["guide_counts"] = _as_count_csr(adata[:, guide_mask].X, "guide count matrix")
    result.uns["guide_ids"] = np.asarray(adata.var_names[guide_mask].astype(str), dtype=str)
    return result


def _read_cellranger_h5(path: str | Path) -> Any:
    try:
        import scanpy as sc
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError("Scanpy is required for Cell Ranger HDF5 input.") from exc
    return _split_feature_types(sc.read_10x_h5(path, gex_only=False))


def _read_h5ad(path: str | Path) -> Any:
    source = ad.read_h5ad(path)
    if "guide_counts" in source.obsm and "guide_ids" in source.uns:
        source.X = _as_count_csr(source.X, "expression count matrix")
        source.obsm["guide_counts"] = _as_count_csr(
            source.obsm["guide_counts"], "guide count matrix"
        )
        return source
    return _split_feature_types(source)


def _read_csv_fixture(input_config: dict[str, Any]) -> Any:
    counts = read_counts_csv(input_config["expression_counts"])
    guides = read_counts_csv(input_config["guide_counts"])
    if counts.index.has_duplicates or counts.columns.has_duplicates:
        raise ValueError("Expression count matrix contains duplicate features or cells")
    if guides.index.has_duplicates or guides.columns.has_duplicates:
        raise ValueError("Guide count matrix contains duplicate features or cells")
    guides = guides.reindex(columns=counts.columns, fill_value=0)
    result = ad.AnnData(
        X=_as_count_csr(counts.T, "expression count matrix"),
        obs=pd.DataFrame(index=counts.columns),
        var=pd.DataFrame(index=counts.index),
    )
    result.obsm["guide_counts"] = _as_count_csr(guides.T, "guide count matrix")
    result.uns["guide_ids"] = np.asarray(guides.index, dtype=str)
    return result


def _align_guide_matrix(matrix: sparse.csr_matrix, guides: list[str], union: list[str]):
    guide_to_column = {guide: index for index, guide in enumerate(union)}
    mapping = np.asarray([guide_to_column[guide] for guide in guides])
    coo = matrix.tocoo()
    return sparse.csr_matrix(
        (coo.data, (coo.row, mapping[coo.col])),
        shape=(matrix.shape[0], len(union)),
    )


def _read_manifest(path: str | Path) -> Any:
    manifest_path = Path(path)
    manifest = pd.read_csv(manifest_path, keep_default_na=False)
    missing = MANIFEST_REQUIRED - set(manifest.columns)
    if missing:
        raise ValueError(f"Input manifest is missing columns: {sorted(missing)}")
    if manifest.empty or manifest["sample_id"].duplicated().any():
        raise ValueError("Input manifest must contain unique sample_id values")

    datasets = []
    guide_matrices = []
    guide_lists = []
    for row in manifest.to_dict("records"):
        sample_id = str(row["sample_id"]).strip()
        input_path = Path(str(row["input_path"]))
        if not input_path.is_absolute():
            input_path = manifest_path.parent / input_path
        input_format = str(row["input_format"]).strip().lower()
        if input_format == "cellranger_h5":
            dataset = _read_cellranger_h5(input_path)
        elif input_format == "h5ad":
            dataset = _read_h5ad(input_path)
        else:
            raise ValueError(f"Unsupported manifest input_format for {sample_id}: {input_format}")

        metadata_path = str(row.get("metadata_path", "")).strip()
        if metadata_path:
            resolved = Path(metadata_path)
            if not resolved.is_absolute():
                resolved = manifest_path.parent / resolved
            metadata = pd.read_csv(resolved)
        else:
            metadata = dataset.obs.reset_index(names="cell_barcode")
        defaults = {
            "sample_id": sample_id,
            "donor": row["donor"],
            "batch": row["batch"],
            "condition": row["condition"],
            "cell_type": row.get("cell_type", "unknown"),
        }
        prepared = _prepare_metadata(metadata, dataset.obs_names.astype(str), defaults)
        namespaced = pd.Index([f"{sample_id}#{barcode}" for barcode in dataset.obs_names])
        dataset.obs = prepared.copy()
        dataset.obs_names = namespaced
        dataset.obs.index.name = "cell_barcode"
        guide_matrices.append(sparse.csr_matrix(dataset.obsm.pop("guide_counts")))
        guide_lists.append([str(value) for value in dataset.uns.pop("guide_ids")])
        datasets.append(dataset)

    guide_union = sorted(set().union(*map(set, guide_lists)))
    aligned_guides = [
        _align_guide_matrix(matrix, guides, guide_union)
        for matrix, guides in zip(guide_matrices, guide_lists, strict=True)
    ]
    combined = ad.concat(datasets, join="outer", merge="same", fill_value=0, index_unique=None)
    combined.X = _as_count_csr(combined.X, "combined expression count matrix")
    combined.obsm["guide_counts"] = sparse.vstack(aligned_guides, format="csr")
    combined.uns["guide_ids"] = np.asarray(guide_union, dtype=str)
    return combined


def load_analysis(path: str | Path) -> Any:
    """Load and validate the canonical sparse analysis object."""
    adata = ad.read_h5ad(path)
    if "guide_counts" not in adata.obsm or "guide_ids" not in adata.uns:
        raise ValueError("Analysis H5AD is missing obsm['guide_counts'] or uns['guide_ids']")
    adata.X = _as_count_csr(adata.X, "expression count matrix")
    adata.obsm["guide_counts"] = _as_count_csr(adata.obsm["guide_counts"], "guide count matrix")
    return adata


def prepare_inputs(
    config: dict[str, Any],
    adata_out: str | Path,
    matrix_out: str | Path,
    genes_out: str | Path,
    barcodes_out: str | Path,
    metadata_out: str | Path,
) -> None:
    """Normalize supported inputs into a sparse, analysis-ready H5AD contract."""
    input_config = config["input"]
    mode = input_config["mode"]
    if mode == "csv":
        adata = _read_csv_fixture(input_config)
    elif mode == "cellranger_h5":
        adata = _read_cellranger_h5(input_config["cellranger_h5"])
    elif mode == "h5ad":
        adata = _read_h5ad(input_config["h5ad"])
    elif mode == "manifest":
        adata = _read_manifest(input_config["manifest"])
    else:
        raise ValueError(f"Unsupported input mode: {mode}")

    metadata_path = input_config.get("metadata")
    if mode != "manifest":
        metadata = (
            pd.read_csv(metadata_path)
            if metadata_path
            else adata.obs.reset_index(names="cell_barcode")
        )
        adata.obs = _prepare_metadata(
            metadata,
            adata.obs_names.astype(str),
            input_config.get("sample_metadata"),
        )
    sample_design = input_config.get("sample_design")
    if sample_design:
        _validate_sample_design(adata.obs, sample_design)
    adata.obs.index.name = "cell_barcode"
    adata.layers["counts"] = adata.X.copy()
    adata.uns["pipeline_schema_version"] = "1.0"

    adata.write_h5ad(ensure_parent(adata_out), compression="gzip")
    mmwrite(ensure_parent(matrix_out), sparse.csc_matrix(adata.X.T), field="integer")
    ensure_parent(genes_out).write_text("\n".join(adata.var_names.astype(str)) + "\n")
    ensure_parent(barcodes_out).write_text("\n".join(adata.obs_names.astype(str)) + "\n")
    write_table(adata.obs.reset_index(), metadata_out)
