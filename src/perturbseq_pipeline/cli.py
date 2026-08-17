import argparse
import json
from pathlib import Path

import pandas as pd

from perturbseq_pipeline.demo import make_demo
from perturbseq_pipeline.design_gate import evaluate_design_gate, require_design_gate
from perturbseq_pipeline.expression_diagnostics import (
    contingency_diagnostics,
    expression_pca_sparse,
    write_recommendations,
)
from perturbseq_pipeline.guide_diagnostics import guide_batch_diagnostics
from perturbseq_pipeline.guides import (
    assign_guides_sparse,
    summarize_guide_qc,
    summarize_low_moi_qc,
)
from perturbseq_pipeline.inputs import (
    as_bool,
    load_analysis,
    prepare_inputs,
    read_guides,
    read_yaml,
    write_matrix,
    write_table,
)
from perturbseq_pipeline.provenance import write_run_manifest
from perturbseq_pipeline.pseudobulk import aggregate_pseudobulk_sparse
from perturbseq_pipeline.qc import (
    calculate_qc_metrics,
    call_cells,
    filter_cells,
    score_doublets_scrublet,
)
from perturbseq_pipeline.reporting import (
    save_batch_plots,
    save_guide_batch_plots,
    save_qc_plots,
    write_html_report,
)


def cmd_make_demo(args: argparse.Namespace) -> None:
    make_demo(args.outdir)


def cmd_prepare_inputs(args: argparse.Namespace) -> None:
    prepare_inputs(
        read_yaml(args.config),
        adata_out=args.adata_out,
        matrix_out=args.matrix_out,
        genes_out=args.genes_out,
        barcodes_out=args.barcodes_out,
        metadata_out=args.metadata_out,
    )


def cmd_qc(args: argparse.Namespace) -> None:
    cfg = read_yaml(args.config)
    adata = load_analysis(args.adata)
    qc_cfg = cfg["qc"]
    qc = calculate_qc_metrics(adata)
    calling_cfg = cfg.get("cell_calling", {"method": "filtered_input"})
    qc = call_cells(
        qc,
        method=calling_cfg.get("method", "filtered_input"),
        min_counts=calling_cfg.get("min_counts", 100),
        min_genes=calling_cfg.get("min_genes", 50),
    )
    qc = filter_cells(
        qc,
        min_counts=qc_cfg["min_counts"],
        min_genes=qc_cfg["min_genes"],
        max_mito_fraction=qc_cfg["max_mito_fraction"],
    )
    qc = score_doublets_scrublet(
        adata,
        qc,
        expected_doublet_rate=qc_cfg["expected_doublet_rate"],
        threshold=qc_cfg.get("doublet_score_threshold"),
    )
    write_table(qc, args.metrics)
    save_qc_plots(qc, args.plots)


def cmd_guides(args: argparse.Namespace) -> None:
    cfg = read_yaml(args.config)
    guide_cfg = cfg["guide_assignment"]
    adata = load_analysis(args.adata)
    guides = read_guides(args.guide_reference)
    assignments = assign_guides_sparse(
        adata.obsm["guide_counts"],
        adata.obs_names.astype(str).tolist(),
        [str(guide) for guide in adata.uns["guide_ids"]],
        guides,
        min_umi=guide_cfg["min_umi"],
        min_fraction=guide_cfg["min_fraction"],
    )
    write_table(assignments, args.assignments)
    write_table(summarize_guide_qc(assignments), args.summary)
    write_table(
        summarize_low_moi_qc(
            assignments,
            expected_moi=cfg["experiment"]["expected_moi"],
            max_observed_multiplet_fraction=guide_cfg["max_observed_multiplet_fraction"],
        ),
        args.low_moi_summary,
    )


def cmd_batch(args: argparse.Namespace) -> None:
    cfg = read_yaml(args.config)
    metadata = pd.read_csv(args.metadata)
    batch_cfg = cfg["batch"]
    eligibility_column = batch_cfg.get("eligibility_column")
    if eligibility_column:
        eligible = as_bool(metadata[eligibility_column])
        metadata = metadata.loc[eligible].copy()
    adata = load_analysis(args.adata)
    pca_scores, variance_ratio = expression_pca_sparse(
        adata,
        metadata,
        highly_variable_genes=cfg["qc"].get("highly_variable_genes", 2000),
    )
    diagnostics = contingency_diagnostics(
        metadata,
        batch_cfg["covariates"],
        pca_scores=pca_scores,
        variance_ratio=variance_ratio,
    )
    write_table(diagnostics, args.metrics)
    Path(args.recommendations).write_text(write_recommendations(diagnostics))
    save_batch_plots(pca_scores, metadata, batch_cfg["covariates"], args.plots)


def cmd_guide_batch(args: argparse.Namespace) -> None:
    cfg = read_yaml(args.config)
    guide_batch_cfg = cfg["guide_batch"]
    assignments = pd.read_csv(args.assignments)
    metadata = pd.read_csv(args.metadata)
    by_sample, tests, composition = guide_batch_diagnostics(
        assignments,
        metadata,
        sample_column=guide_batch_cfg.get("sample_column", "sample_id"),
        batch_column=guide_batch_cfg.get("batch_column", "batch"),
        min_samples_per_batch=guide_batch_cfg.get("min_samples_per_batch", 2),
        fdr_threshold=guide_batch_cfg.get("fdr_threshold", 0.05),
        effect_size_threshold=guide_batch_cfg.get("effect_size_threshold", 0.10),
    )
    write_table(by_sample, args.by_sample)
    write_table(tests, args.tests)
    write_table(composition, args.composition)
    save_guide_batch_plots(
        by_sample,
        composition,
        guide_batch_cfg.get("batch_column", "batch"),
        args.plots,
    )


def cmd_design_gate(args: argparse.Namespace) -> None:
    cfg = read_yaml(args.config)
    metadata = pd.read_csv(args.metadata)
    mixscape = pd.read_csv(args.mixscape)[
        ["cell_barcode", "mixscape_class", "escaped_perturbation"]
    ]
    metadata = metadata.merge(mixscape, on="cell_barcode", how="left", validate="one_to_one")
    gate_cfg = cfg["design_gate"]
    diagnostics, decision = evaluate_design_gate(
        metadata,
        batch_column=gate_cfg.get("batch_column", "batch"),
        confounders=gate_cfg.get(
            "confounders",
            ["donor", "assigned_guide", "target_gene", "escaped_perturbation"],
        ),
        eligibility_column=gate_cfg.get("eligibility_column", "include_primary_analysis"),
    )
    write_table(diagnostics, args.diagnostics)
    output = Path(args.decision)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(decision, indent=2) + "\n")


def cmd_pseudobulk(args: argparse.Namespace) -> None:
    cfg = read_yaml(args.config)
    if args.design_gate:
        require_design_gate(json.loads(Path(args.design_gate).read_text()))
    adata = load_analysis(args.adata)
    metadata = pd.read_csv(args.metadata)
    pseudobulk_config = dict(cfg["pseudobulk"])
    if args.mixscape_scores:
        mixscape = pd.read_csv(args.mixscape_scores)[["cell_barcode", "mixscape_class"]]
        metadata = metadata.merge(mixscape, on="cell_barcode", how="left", validate="one_to_one")
        metadata["include_ko_sensitivity"] = as_bool(
            metadata["include_primary_analysis"]
        ) & metadata["mixscape_class"].isin(["KO", "NT"])
        pseudobulk_config["eligibility_column"] = "include_ko_sensitivity"
    pseudobulk, design = aggregate_pseudobulk_sparse(adata, metadata, **pseudobulk_config)
    write_matrix(pseudobulk, args.counts_out)
    write_table(design, args.design_out)


def cmd_report(args: argparse.Namespace) -> None:
    write_html_report(
        output_path=args.output,
        qc_path=args.qc,
        low_moi_path=args.low_moi,
        batch_path=args.batch,
        guide_batch_sample_path=args.guide_batch_sample,
        guide_batch_tests_path=args.guide_batch_tests,
        guide_batch_composition_path=args.guide_batch_composition,
        gate_diagnostics_path=args.gate_diagnostics,
        gate_decision_path=args.gate_decision,
        mixscape_path=args.mixscape,
        contrast_path=args.contrasts,
        de_path=args.de,
        sensitivity_contrast_path=args.sensitivity_contrasts,
        sensitivity_de_path=args.sensitivity_de,
        pathway_path=args.pathway,
        top_n=read_yaml(args.config)["reports"]["top_n_genes"],
    )


def cmd_provenance(args: argparse.Namespace) -> None:
    write_run_manifest(args.config, args.adata, args.output)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="perturbseq-pipeline")
    sub = parser.add_subparsers(required=True)

    p = sub.add_parser("make-demo")
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=cmd_make_demo)

    p = sub.add_parser("prepare-inputs")
    p.add_argument("--config", required=True)
    p.add_argument("--adata-out", required=True)
    p.add_argument("--matrix-out", required=True)
    p.add_argument("--genes-out", required=True)
    p.add_argument("--barcodes-out", required=True)
    p.add_argument("--metadata-out", required=True)
    p.set_defaults(func=cmd_prepare_inputs)

    p = sub.add_parser("qc")
    p.add_argument("--config", required=True)
    p.add_argument("--adata", required=True)
    p.add_argument("--metrics", required=True)
    p.add_argument("--plots", required=True)
    p.set_defaults(func=cmd_qc)

    p = sub.add_parser("guides")
    p.add_argument("--config", required=True)
    p.add_argument("--adata", required=True)
    p.add_argument("--guide-reference", required=True)
    p.add_argument("--assignments", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--low-moi-summary", required=True)
    p.set_defaults(func=cmd_guides)

    p = sub.add_parser("batch")
    p.add_argument("--config", required=True)
    p.add_argument("--adata", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--metrics", required=True)
    p.add_argument("--recommendations", required=True)
    p.add_argument("--plots", required=True)
    p.set_defaults(func=cmd_batch)

    p = sub.add_parser("guide-batch")
    p.add_argument("--config", required=True)
    p.add_argument("--assignments", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--by-sample", required=True)
    p.add_argument("--tests", required=True)
    p.add_argument("--composition", required=True)
    p.add_argument("--plots", required=True)
    p.set_defaults(func=cmd_guide_batch)

    p = sub.add_parser("design-gate")
    p.add_argument("--config", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--mixscape", required=True)
    p.add_argument("--diagnostics", required=True)
    p.add_argument("--decision", required=True)
    p.set_defaults(func=cmd_design_gate)

    p = sub.add_parser("pseudobulk")
    p.add_argument("--config", required=True)
    p.add_argument("--adata", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--mixscape-scores")
    p.add_argument("--design-gate")
    p.add_argument("--counts-out", required=True)
    p.add_argument("--design-out", required=True)
    p.set_defaults(func=cmd_pseudobulk)

    p = sub.add_parser("report")
    p.add_argument("--config", required=True)
    p.add_argument("--qc", required=True)
    p.add_argument("--low-moi", required=True)
    p.add_argument("--batch", required=True)
    p.add_argument("--guide-batch-sample", required=True)
    p.add_argument("--guide-batch-tests", required=True)
    p.add_argument("--guide-batch-composition", required=True)
    p.add_argument("--gate-diagnostics", required=True)
    p.add_argument("--gate-decision", required=True)
    p.add_argument("--mixscape", required=True)
    p.add_argument("--contrasts", required=True)
    p.add_argument("--de", required=True)
    p.add_argument("--sensitivity-contrasts", required=True)
    p.add_argument("--sensitivity-de", required=True)
    p.add_argument("--pathway", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_report)

    p = sub.add_parser("provenance")
    p.add_argument("--config", required=True)
    p.add_argument("--adata", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_provenance)

    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
