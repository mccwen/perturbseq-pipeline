import html
import json
from pathlib import Path

import numpy as np
import pandas as pd


def save_qc_plots(qc: pd.DataFrame, output_pdf: str | Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(output_pdf)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
    axes[0].hist(qc["total_counts"], bins=30, color="#2f6f8f")
    axes[1].hist(qc["detected_genes"], bins=30, color="#3f7d5c")
    axes[2].hist(qc["mito_fraction"], bins=30, color="#b45f3c")
    axes[0].set_title("Library size")
    axes[1].set_title("Detected genes")
    axes[2].set_title("Mito fraction")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def save_batch_plots(
    pca_scores: pd.DataFrame,
    metadata: pd.DataFrame,
    covariates: list[str],
    output_pdf: str | Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    plot_data = pca_scores.merge(metadata, on="cell_barcode", validate="one_to_one")
    out = Path(output_pdf)
    out.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out) as pdf:
        for covariate in covariates:
            if covariate not in plot_data:
                continue
            fig, ax = plt.subplots(figsize=(6.2, 5.2))
            levels = sorted(plot_data[covariate].astype(str).unique())
            palette = plt.get_cmap("tab10")
            for index, level in enumerate(levels):
                subset = plot_data[covariate].astype(str).eq(level)
                ax.scatter(
                    plot_data.loc[subset, "PC1"],
                    plot_data.loc[subset, "PC2"],
                    s=12,
                    alpha=0.65,
                    color=palette(index % 10),
                    label=level,
                )
            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
            ax.set_title(f"Expression PCA by {covariate}")
            ax.legend(frameon=False, fontsize=8)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def save_guide_batch_plots(
    by_sample: pd.DataFrame,
    composition: pd.DataFrame,
    batch_column: str,
    output_pdf: str | Path,
) -> None:
    """Plot replicate-level guide QC and batch-stratified guide composition."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    metrics = [
        "guide_positive_rate",
        "singlet_rate",
        "unassigned_rate",
        "multiplet_rate_among_guide_positive",
        "mean_guide_umi_total",
        "median_guide_umi_positive",
        "median_top_guide_fraction",
    ]
    batches = sorted(by_sample[batch_column].astype(str).unique())
    palette = plt.get_cmap("tab10")
    out = Path(output_pdf)
    out.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out) as pdf:
        fig, axes = plt.subplots(3, 3, figsize=(11, 10))
        for ax, metric in zip(axes.flat, metrics):
            groups = [
                by_sample.loc[by_sample[batch_column].astype(str).eq(batch), metric]
                .dropna()
                .to_numpy()
                for batch in batches
            ]
            ax.boxplot(groups, tick_labels=batches, widths=0.55, showfliers=False)
            for index, values in enumerate(groups, start=1):
                if len(values) == 0:
                    continue
                offsets = np.linspace(-0.08, 0.08, len(values))
                ax.scatter(
                    index + offsets,
                    values,
                    s=24,
                    color=palette((index - 1) % 10),
                    alpha=0.8,
                    zorder=3,
                )
            ax.set_title(metric.replace("_", " "))
            ax.set_xlabel(batch_column)
        for ax in axes.flat[len(metrics) :]:
            ax.set_visible(False)
        fig.suptitle("Guide-capture QC by biological sample", fontsize=13)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        pivot = composition.pivot(
            index=batch_column,
            columns="assigned_guide",
            values="guide_fraction_among_singlets",
        ).fillna(0)
        fig, ax = plt.subplots(figsize=(8, 5.2))
        pivot.plot(kind="bar", stacked=True, ax=ax, colormap="tab20", width=0.7)
        ax.set_ylabel("Fraction among guide singlets")
        ax.set_xlabel(batch_column)
        ax.set_title("Assigned-guide composition by batch")
        ax.legend(title="Guide", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def _table_html(frame: pd.DataFrame, columns: list[str], limit: int = 20) -> str:
    available = [column for column in columns if column in frame]
    if not available or frame.empty:
        return "<p>No results available.</p>"
    display = frame.loc[:, available].head(limit).copy()
    table = display.to_html(
        index=False,
        border=0,
        escape=True,
        na_rep="",
        float_format=lambda value: f"{value:.4g}",
    )
    return f"<div class='table-wrap'>{table}</div>"


def write_html_report(
    output_path: str | Path,
    qc_path: str | Path,
    low_moi_path: str | Path,
    batch_path: str | Path,
    guide_batch_sample_path: str | Path,
    guide_batch_tests_path: str | Path,
    guide_batch_composition_path: str | Path,
    gate_diagnostics_path: str | Path,
    gate_decision_path: str | Path,
    mixscape_path: str | Path,
    contrast_path: str | Path,
    de_path: str | Path,
    sensitivity_contrast_path: str | Path,
    sensitivity_de_path: str | Path,
    pathway_path: str | Path,
    top_n: int = 25,
) -> None:
    qc = pd.read_csv(qc_path)
    low_moi = pd.read_csv(low_moi_path)
    batch = pd.read_csv(batch_path)
    guide_batch_sample = pd.read_csv(guide_batch_sample_path)
    guide_batch_tests = pd.read_csv(guide_batch_tests_path)
    guide_batch_composition = pd.read_csv(guide_batch_composition_path)
    gate_diagnostics = pd.read_csv(gate_diagnostics_path)
    gate_decision = json.loads(Path(gate_decision_path).read_text())
    mixscape = pd.read_csv(mixscape_path)
    contrasts = pd.read_csv(contrast_path)
    de = pd.read_csv(de_path)
    sensitivity_contrasts = pd.read_csv(sensitivity_contrast_path)
    sensitivity_de = pd.read_csv(sensitivity_de_path)
    pathway = pd.read_csv(pathway_path)

    low = low_moi.iloc[0]
    tested = contrasts[contrasts["status"] == "tested"]
    significant = de[(de["FDR"] < 0.05) & (de["logFC"].abs() >= 0.5)]
    top_de = de.sort_values(["FDR", "PValue"]).head(top_n)
    top_pathways = pathway.sort_values(["padj", "pval"]).head(top_n)
    class_counts = (
        mixscape["mixscape_class"].value_counts().rename_axis("class").reset_index(name="cells")
    )

    metrics = [
        ("Cells", len(qc)),
        ("QC pass", int(qc["pass_qc"].sum())),
        ("Guide singlets", int(low["n_single_guide"])),
        ("Guide multiplets", int(low["n_guide_multiplet"])),
        ("Guide batch flags", int(guide_batch_tests["batch_effect_detected"].sum())),
        ("Design gate", gate_decision["status"].upper()),
        ("Tested contrasts", len(tested)),
        ("Significant genes", len(significant)),
        ("KO-only contrasts", int(sensitivity_contrasts["status"].eq("tested").sum())),
    ]
    metric_html = "".join(
        f"<div class='metric'><span>{html.escape(label)}</span><strong>{value}</strong></div>"
        for label, value in metrics
    )

    sections = [
        (
            "Low-MOI QC",
            _table_html(
                low_moi,
                [
                    "expected_moi",
                    "apparent_moi_from_guide_detection",
                    "observed_guide_positive_fraction",
                    "observed_multiplet_fraction_among_guide_positive",
                    "poisson_expected_multiplet_fraction_at_target_moi",
                    "passes_low_moi_qc",
                ],
                1,
            ),
        ),
        (
            "Expression Batch Diagnostics",
            _table_html(
                batch,
                [
                    "covariate",
                    "n_levels",
                    "zero_fraction",
                    "weighted_pc_r2",
                    "max_pc_r2",
                    "recommendation",
                ],
            ),
        ),
        (
            "Guide-Capture Batch Tests",
            _table_html(
                guide_batch_tests.sort_values(
                    ["batch_effect_detected", "fdr"], ascending=[False, True]
                ),
                [
                    "scope",
                    "metric",
                    "n_batches",
                    "n_samples",
                    "min_samples_per_batch",
                    "effect_size_epsilon_squared",
                    "p_value",
                    "fdr",
                    "status",
                    "batch_effect_detected",
                ],
            ),
        ),
        (
            "Guide QC by Sample",
            _table_html(
                guide_batch_sample,
                [
                    "sample_id",
                    "batch",
                    "n_cells",
                    "guide_positive_rate",
                    "singlet_rate",
                    "unassigned_rate",
                    "multiplet_rate_among_guide_positive",
                    "mean_guide_umi_total",
                    "median_guide_umi_positive",
                    "median_top_guide_fraction",
                ],
            ),
        ),
        (
            "Guide Composition by Batch",
            _table_html(
                guide_batch_composition,
                [
                    "batch",
                    "assigned_guide",
                    "n_assigned_cells",
                    "n_single_cells",
                    "guide_fraction_among_singlets",
                ],
            ),
        ),
        (
            "Identifiability Gate",
            f"<p><strong>{html.escape(gate_decision['status'].upper())}:</strong> "
            f"{html.escape(gate_decision['reason'])}</p>"
            + _table_html(
                gate_diagnostics,
                [
                    "batch_column",
                    "confounder",
                    "n_batch_levels",
                    "n_confounder_levels",
                    "connected_components",
                    "completely_confounded",
                    "status",
                ],
            ),
        ),
        ("Mixscape Classes", _table_html(class_counts, ["class", "cells"])),
        (
            "Differential Expression Contrasts",
            _table_html(
                contrasts,
                [
                    "stratum",
                    "target_gene",
                    "reference",
                    "n_target_samples",
                    "n_reference_samples",
                    "status",
                    "reason",
                    "model_formula",
                ],
            ),
        ),
        (
            "Top Differentially Expressed Genes",
            _table_html(
                top_de,
                ["gene", "target_gene", "reference", "stratum", "logFC", "F", "PValue", "FDR"],
                top_n,
            ),
        ),
        (
            "KO-Only Sensitivity Analysis",
            _table_html(
                sensitivity_de.sort_values(["FDR", "PValue"]),
                ["gene", "target_gene", "reference", "stratum", "logFC", "F", "PValue", "FDR"],
                top_n,
            ),
        ),
        (
            "Top Pathways",
            _table_html(
                top_pathways,
                ["pathway", "target_gene", "reference", "stratum", "NES", "pval", "padj", "size"],
                top_n,
            ),
        ),
    ]
    section_html = "".join(
        f"<section><h2>{title}</h2>{content}</section>" for title, content in sections
    )

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Low-MOI Perturb-seq Summary</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 0; color: #17212b; background: #f5f7f9; }}
    header {{ padding: 2rem max(1.25rem, 5vw); color: white; background: #17324d; }}
    h1 {{ margin: 0; font-size: 1.8rem; }}
    header p {{ margin: .45rem 0 0; color: #dbe7f2; }}
    main {{ max-width: 1200px; margin: 0 auto; padding: 1.5rem; }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: .75rem; }}
    .metric {{ padding: 1rem; background: white; border: 1px solid #d9e0e6; border-radius: 6px; }}
    .metric span {{ display: block; color: #536273; font-size: .8rem; }}
    .metric strong {{ display: block; margin-top: .3rem; font-size: 1.5rem; }}
    section {{ margin-top: 1.5rem; padding: 1rem; background: white; border-top: 3px solid #2f6f8f; }}
    h2 {{ margin: 0 0 .8rem; font-size: 1.1rem; }}
    .table-wrap {{ overflow-x: auto; }}
    table {{ width: 100%; border-collapse: collapse; font-size: .82rem; }}
    th, td {{ padding: .55rem; border-bottom: 1px solid #e3e8ed; text-align: left; white-space: nowrap; }}
    th {{ color: #324456; background: #eef3f6; }}
  </style>
</head>
<body>
  <header><h1>Low-MOI Perturb-seq Analysis</h1><p>QC, guide assignment, perturbation efficacy, differential expression, and pathways</p></header>
  <main><div class="metrics">{metric_html}</div>{section_html}</main>
</body>
</html>
"""
    out.write_text(document, encoding="utf-8")
