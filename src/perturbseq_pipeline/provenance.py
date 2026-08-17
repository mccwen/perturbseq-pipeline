import hashlib
import json
import subprocess
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from perturbseq_pipeline.inputs import (
    configured_input_paths,
    ensure_parent,
    load_analysis,
    read_yaml,
    resolve_path,
)

TRACKED_PACKAGES = [
    "perturbseq-pipeline",
    "anndata",
    "matplotlib",
    "numpy",
    "pandas",
    "scanpy",
    "scipy",
    "scrublet",
    "snakemake",
]


def sha256_file(path: str | Path) -> str:
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _portable_path(path: str | Path, root: str | Path) -> str:
    """Represent repository files without exposing a local absolute path."""
    source = Path(path).resolve()
    try:
        relative = source.relative_to(Path(root).resolve())
    except ValueError:
        return f"<external>/{source.name}"
    return f"./{relative.as_posix()}"


def _portable_config(value: Any, root: Path) -> Any:
    """Replace absolute paths nested in a resolved configuration."""
    if isinstance(value, dict):
        return {key: _portable_config(item, root) for key, item in value.items()}
    if isinstance(value, list):
        return [_portable_config(item, root) for item in value]
    if isinstance(value, str) and Path(value).is_absolute():
        return _portable_path(value, root)
    return value


def file_fingerprint(path: str | Path, root: str | Path = ".") -> dict[str, Any]:
    source = Path(path).resolve()
    stat = source.stat()
    return {
        "path": _portable_path(source, root),
        "size_bytes": stat.st_size,
        "sha256": sha256_file(source),
    }


def configuration_root(config: dict[str, Any], config_path: str | Path) -> Path:
    """Locate the project root for source or workflow-resolved configurations."""
    config_path = Path(config_path).resolve()
    input_config = config.get("input_config")
    if input_config:
        source = Path(input_config)
        if source.is_absolute():
            return source.parent.parent
        for candidate in [Path.cwd().resolve(), *config_path.parents]:
            if (candidate / source).is_file():
                return candidate
    return config_path.parent.parent


def configured_source_paths(config: dict[str, Any], config_path: str | Path) -> list[Path]:
    """Resolve all configured primary inputs and biological reference files."""
    root = configuration_root(config, config_path)
    sources = configured_input_paths(config, root)
    sources.extend(
        [
            resolve_path(config["guide_assignment"]["reference"], root),
            resolve_path(config["pathway_analysis"]["gene_sets"], root),
        ]
    )
    return list(dict.fromkeys(sources))


def _git_value(arguments: list[str], root: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", *arguments],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def _git_state(root: Path) -> dict[str, Any]:
    commit = _git_value(["rev-parse", "HEAD"], root)
    status = _git_value(["status", "--porcelain"], root)
    return {
        "commit": commit,
        "branch": _git_value(["branch", "--show-current"], root),
        "dirty": None if status is None else bool(status),
    }


def _package_versions() -> dict[str, str | None]:
    versions = {}
    for package in TRACKED_PACKAGES:
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = None
    return versions


def write_run_manifest(
    config_path: str | Path,
    adata_path: str | Path,
    output_path: str | Path,
) -> dict[str, Any]:
    config_path = Path(config_path).resolve()
    config = read_yaml(config_path)
    root = configuration_root(config, config_path)
    adata = load_analysis(adata_path)
    guide_matrix = adata.obsm["guide_counts"]
    manifest = {
        "schema_version": "1.1",
        "workflow": "post-cell-ranger-low-moi-perturb-seq",
        "source_control": _git_state(root),
        "software": _package_versions(),
        "configuration": {
            "path": _portable_path(config_path, root),
            "sha256": sha256_file(config_path),
            "resolved": _portable_config(config, root),
        },
        "source_inputs": [
            file_fingerprint(path, root) for path in configured_source_paths(config, config_path)
        ],
        "analysis_object": {
            "path": _portable_path(adata_path, root),
            "sha256": sha256_file(adata_path),
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "n_guides": int(guide_matrix.shape[1]),
            "expression_nonzero": int(adata.X.nnz),
            "guide_nonzero": int(guide_matrix.nnz),
            "pipeline_schema_version": adata.uns.get("pipeline_schema_version"),
        },
    }
    output = ensure_parent(output_path)
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest
