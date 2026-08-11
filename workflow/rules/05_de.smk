rule differential_expression:
    input:
        counts="data/processed/04_pseudobulk/pseudobulk_counts.csv",
        design="data/processed/04_pseudobulk/design_matrix.csv"
    output:
        results="data/processed/05_de/de_results.csv",
        plots="data/processed/05_de/de_plots.pdf"
    script:
        "../../scripts/run_edger.R"

