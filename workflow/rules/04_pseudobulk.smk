rule pseudobulk:
    input:
        counts="data/demo/expression_counts.csv",
        metadata="data/processed/03_batch/cell_metadata_with_guides.csv",
        config="config/config.yaml"
    output:
        counts_out="data/processed/04_pseudobulk/pseudobulk_counts.csv",
        design_out="data/processed/04_pseudobulk/design_matrix.csv"
    shell:
        "python -m perturbseq_pipeline pseudobulk "
        "--config {input.config} "
        "--counts {input.counts} "
        "--metadata {input.metadata} "
        "--counts-out {output.counts_out} "
        "--design-out {output.design_out}"

