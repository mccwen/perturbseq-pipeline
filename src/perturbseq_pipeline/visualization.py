from __future__ import annotations

from pathlib import Path

import pandas as pd


def save_qc_plots(qc: pd.DataFrame, output_pdf: str | Path) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    out = Path(output_pdf)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
    sns.histplot(qc["total_counts"], ax=axes[0], bins=30)
    sns.histplot(qc["detected_genes"], ax=axes[1], bins=30)
    sns.histplot(qc["mito_fraction"], ax=axes[2], bins=30)
    axes[0].set_title("Library size")
    axes[1].set_title("Detected genes")
    axes[2].set_title("Mito fraction")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def write_html_report(output_path: str | Path, sections: dict[str, str]) -> None:
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    body = "\n".join(f"<h2>{title}</h2>\n<p>{text}</p>" for title, text in sections.items())
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Perturb-seq Summary</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; line-height: 1.5; color: #1f2933; }}
    h1, h2 {{ color: #12355b; }}
    code {{ background: #eef2f7; padding: 0.15rem 0.3rem; border-radius: 4px; }}
  </style>
</head>
<body>
  <h1>Perturb-seq / CROP-seq Analysis Summary</h1>
  {body}
</body>
</html>
"""
    out.write_text(html)
