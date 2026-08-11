from __future__ import annotations

import numpy as np
import pandas as pd


def classify_escaped_perturbations(
    assignments: pd.DataFrame,
    expression_signature: pd.Series | None = None,
    threshold: float = 0.70,
) -> pd.DataFrame:
    """Create a Mixscape-like KO / NP / NT classification table.

    The R script provides the production Mixscape hook. This Python fallback makes the
    workflow testable and documents the expected output schema.
    """
    scores = assignments[["cell_barcode", "target_gene", "assigned_guide", "is_non_targeting"]].copy()
    if expression_signature is not None:
        scores = scores.join(expression_signature.rename("perturbation_score"), on="cell_barcode")
    else:
        raw = assignments["guide_fraction"].fillna(0).to_numpy()
        scores["perturbation_score"] = np.clip(raw, 0, 1)

    scores["mixscape_class"] = np.select(
        [
            scores["is_non_targeting"],
            scores["perturbation_score"] >= threshold,
        ],
        ["NT", "KO"],
        default="NP",
    )
    scores["escaped_perturbation"] = scores["mixscape_class"].eq("NP")
    return scores

