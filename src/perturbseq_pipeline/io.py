from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


def read_yaml(path: str | Path) -> dict[str, Any]:
    try:
        import yaml
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "PyYAML is required to read config files. Install the project environment with "
            "`conda env create -f environment.yml` or `pip install pyyaml`."
        ) from exc
    with Path(path).open() as handle:
        return yaml.safe_load(handle)


def ensure_parent(path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    return out


def read_counts_csv(path: str | Path) -> pd.DataFrame:
    """Read a genes x cells count matrix stored as CSV."""
    counts = pd.read_csv(path, index_col=0)
    counts.index = counts.index.astype(str)
    counts.columns = counts.columns.astype(str)
    return counts


def write_table(df: pd.DataFrame, path: str | Path) -> None:
    out = ensure_parent(path)
    df.to_csv(out, index=False)


def write_matrix(df: pd.DataFrame, path: str | Path) -> None:
    out = ensure_parent(path)
    df.to_csv(out)


def read_samples(path: str | Path) -> pd.DataFrame:
    samples = pd.read_csv(path)
    required = {"sample_id", "donor", "batch", "condition", "cellranger_matrix", "guide_matrix"}
    missing = required - set(samples.columns)
    if missing:
        raise ValueError(f"samples file is missing columns: {sorted(missing)}")
    return samples


def read_guides(path: str | Path) -> pd.DataFrame:
    guides = pd.read_csv(path)
    required = {"guide_id", "target_gene", "is_non_targeting"}
    missing = required - set(guides.columns)
    if missing:
        raise ValueError(f"guide reference is missing columns: {sorted(missing)}")
    guides["is_non_targeting"] = guides["is_non_targeting"].astype(str).str.lower().isin(["true", "1", "yes"])
    return guides
