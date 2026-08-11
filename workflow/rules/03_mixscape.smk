rule mixscape_diagnosis:
    input:
        assignments="data/processed/02_guides/guide_assignments.csv",
        config="config/config.yaml"
    output:
        scores="data/processed/03_mixscape/mixscape_scores.csv"
    shell:
        "python -m perturbseq_pipeline mixscape "
        "--config {input.config} "
        "--assignments {input.assignments} "
        "--scores {output.scores}"

