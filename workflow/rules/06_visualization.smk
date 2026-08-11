rule pathway_visualization:
    input:
        de="data/processed/05_de/de_results.csv"
    output:
        pathway="data/processed/06_visualization/pathway_results.csv",
        plots="data/processed/06_visualization/pathway_plots.pdf"
    script:
        "../../scripts/run_fgsea_placeholder.R"

rule report:
    input:
        qc="data/processed/01_qc/qc_metrics.csv",
        de="data/processed/05_de/de_results.csv",
        pathway="data/processed/06_visualization/pathway_results.csv"
    output:
        html="results/final_reports/summary_dashboard.html"
    shell:
        "python -m perturbseq_pipeline report --output {output.html}"

