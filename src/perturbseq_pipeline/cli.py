from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from perturbseq_pipeline.batch_diagnosis import contingency_diagnostics, write_recommendations
from perturbseq_pipeline.demo import make_demo
from perturbseq_pipeline.guide_calling import assign_guides, summarize_guide_qc
from perturbseq_pipeline.io import read_counts_csv, read_guides, read_yaml, write_matrix, write_table
from perturbseq_pipeline.mixscape import classify_escaped_perturbations
from perturbseq_pipeline.pseudobulk import aggregate_pseudobulk
from perturbseq_pipeline.qc import calculate_qc_metrics, filter_cells, flag_doublets
from perturbseq_pipeline.visualization import save_qc_plots, write_html_report


def cmd_make_demo(args: argparse.Namespace) -> None:
    make_demo(args.outdir)


def cmd_qc(args: argparse.Namespace) -> None:
    cfg = read_yaml(args.config)
    counts = read_counts_csv(args.counts)
    qc_cfg = cfg["qc"]
    qc = calculate_qc_metrics(counts)
    qc = filter_cells(
        qc,
        min_counts=qc_cfg["min_counts"],
        min_genes=qc_cfg["min_genes"],
        max_mito_fraction=qc_cfg["max_mito_fraction"],
    )
    qc = flag_doublets(qc, qc_cfg["doublet_score_threshold"])
    write_table(qc, args.metrics)
    save_qc_plots(qc, args.plots)


def cmd_guides(args: argparse.Namespace) -> None:
    cfg = read_yaml(args.config)
    guide_cfg = cfg["guide_assignment"]
    guide_counts = read_counts_csv(args.guide_counts)
    guides = read_guides(args.guide_reference)
    assignments = assign_guides(
        guide_counts,
        guides,
        min_umi=guide_cfg["min_umi"],
        min_fraction=guide_cfg["min_fraction"],
        max_guides_per_cell=guide_cfg["max_guides_per_cell"],
    )
    write_table(assignments, args.assignments)
    write_table(summarize_guide_qc(assignments), args.summary)


def cmd_batch(args: argparse.Namespace) -> None:
    cfg = read_yaml(args.config)
    metadata = pd.read_csv(args.metadata)
    diagnostics = contingency_diagnostics(metadata, cfg["batch"]["covariates"])
    write_table(diagnostics, args.metrics)
    Path(args.recommendations).write_text(write_recommendations(diagnostics))


def cmd_mixscape(args: argparse.Namespace) -> None:
    cfg = read_yaml(args.config)
    assignments = pd.read_csv(args.assignments)
    scores = classify_escaped_perturbations(
        assignments,
        threshold=cfg["mixscape"]["perturbed_probability_threshold"],
    )
    write_table(scores, args.scores)


def cmd_pseudobulk(args: argparse.Namespace) -> None:
    cfg = read_yaml(args.config)
    counts = read_counts_csv(args.counts)
    metadata = pd.read_csv(args.metadata)
    pseudobulk, design = aggregate_pseudobulk(counts, metadata, **cfg["pseudobulk"])
    write_matrix(pseudobulk, args.counts_out)
    write_table(design, args.design_out)


def cmd_report(args: argparse.Namespace) -> None:
    write_html_report(
        args.output,
        {
            "QC": "Cell-level QC metrics and guide summaries are stored under data/processed.",
            "Batch": "Batch recommendations preserve covariates for statistical inference.",
            "DE": "Pseudobulk matrices and design files are ready for edgeR quasi-likelihood testing.",
        },
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="perturbseq-pipeline")
    sub = parser.add_subparsers(required=True)

    p = sub.add_parser("make-demo")
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=cmd_make_demo)

    p = sub.add_parser("qc")
    p.add_argument("--config", required=True)
    p.add_argument("--counts", required=True)
    p.add_argument("--metrics", required=True)
    p.add_argument("--plots", required=True)
    p.set_defaults(func=cmd_qc)

    p = sub.add_parser("guides")
    p.add_argument("--config", required=True)
    p.add_argument("--guide-counts", required=True)
    p.add_argument("--guide-reference", required=True)
    p.add_argument("--assignments", required=True)
    p.add_argument("--summary", required=True)
    p.set_defaults(func=cmd_guides)

    p = sub.add_parser("batch")
    p.add_argument("--config", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--metrics", required=True)
    p.add_argument("--recommendations", required=True)
    p.set_defaults(func=cmd_batch)

    p = sub.add_parser("mixscape")
    p.add_argument("--config", required=True)
    p.add_argument("--assignments", required=True)
    p.add_argument("--scores", required=True)
    p.set_defaults(func=cmd_mixscape)

    p = sub.add_parser("pseudobulk")
    p.add_argument("--config", required=True)
    p.add_argument("--counts", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--counts-out", required=True)
    p.add_argument("--design-out", required=True)
    p.set_defaults(func=cmd_pseudobulk)

    p = sub.add_parser("report")
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_report)

    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
