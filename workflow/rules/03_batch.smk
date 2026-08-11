rule merge_metadata:
    input:
        metadata="data/demo/cell_metadata.csv",
        assignments="data/processed/02_guides/guide_assignments.csv"
    output:
        merged="data/processed/03_batch/cell_metadata_with_guides.csv"
    run:
        import pandas as pd
        meta = pd.read_csv(input.metadata)
        guides = pd.read_csv(input.assignments)
        meta.merge(guides, on="cell_barcode", how="left").to_csv(output.merged, index=False)

rule batch_diagnosis:
    input:
        metadata="data/processed/03_batch/cell_metadata_with_guides.csv",
        config="config/config.yaml"
    output:
        metrics="data/processed/03_batch/batch_metrics.csv",
        recommendations="data/processed/03_batch/batch_recommendations.md"
    shell:
        "python -m perturbseq_pipeline batch "
        "--config {input.config} "
        "--metadata {input.metadata} "
        "--metrics {output.metrics} "
        "--recommendations {output.recommendations}"

