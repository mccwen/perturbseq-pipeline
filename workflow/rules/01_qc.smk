rule make_demo:
    output:
        expression="data/demo/expression_counts.csv",
        guides="data/demo/guide_counts.csv",
        metadata="data/demo/cell_metadata.csv"
    shell:
        "python -m perturbseq_pipeline make-demo --outdir data/demo"

rule qc:
    input:
        counts="data/demo/expression_counts.csv",
        config="config/config.yaml"
    output:
        metrics="data/processed/01_qc/qc_metrics.csv",
        plots="data/processed/01_qc/qc_plots.pdf"
    shell:
        "python -m perturbseq_pipeline qc "
        "--config {input.config} "
        "--counts {input.counts} "
        "--metrics {output.metrics} "
        "--plots {output.plots}"

