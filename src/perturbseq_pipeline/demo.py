from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def make_demo(outdir: str | Path, n_cells: int = 160, seed: int = 7) -> None:
    rng = np.random.default_rng(seed)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    genes = ["TP53", "KRAS", "MYC", "EPCAM", "KRT8", "MKI67", "PTPRC", "MT-CO1", "RPLP0", "ACTB"]
    cells = [f"cell_{i:03d}" for i in range(n_cells)]
    guides = ["NT_001", "NT_002", "TP53_g1", "KRAS_g1", "MYC_g1"]

    expression = pd.DataFrame(rng.negative_binomial(5, 0.35, size=(len(genes), n_cells)), index=genes, columns=cells)
    assigned = rng.choice(guides, size=n_cells, p=[0.18, 0.17, 0.25, 0.20, 0.20])
    guide_counts = pd.DataFrame(0, index=guides, columns=cells)
    for cell, guide in zip(cells, assigned):
        guide_counts.loc[guide, cell] = int(rng.integers(3, 25))

    metadata = pd.DataFrame(
        {
            "cell_barcode": cells,
            "donor": rng.choice(["D01", "D02", "D03", "D04"], size=n_cells),
            "batch": rng.choice(["B1", "B2"], size=n_cells),
            "condition": rng.choice(["vehicle", "stimulated"], size=n_cells),
            "cell_type": rng.choice(["epithelial", "cycling"], size=n_cells),
        }
    )

    expression.to_csv(out / "expression_counts.csv")
    guide_counts.to_csv(out / "guide_counts.csv")
    metadata.to_csv(out / "cell_metadata.csv", index=False)

