rule guide_calling:
    input:
        guide_counts="data/demo/guide_counts.csv",
        guide_reference="config/guides_reference.csv",
        config="config/config.yaml"
    output:
        assignments="data/processed/02_guides/guide_assignments.csv",
        summary="data/processed/02_guides/guide_qc_summary.csv"
    shell:
        "python -m perturbseq_pipeline guides "
        "--config {input.config} "
        "--guide-counts {input.guide_counts} "
        "--guide-reference {input.guide_reference} "
        "--assignments {output.assignments} "
        "--summary {output.summary}"

