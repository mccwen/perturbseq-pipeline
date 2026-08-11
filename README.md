# perturbseq-pipeline

Reproducible Perturb-seq / CROP-seq analysis pipeline from Cell Ranger outputs to biological insight.

This repository is arranged as an interview-ready example for a Computational Scientist role focused on high-content perturbation screens, foundation-model-ready perturbation datasets, and collaborative pipeline development.

## What This Pipeline Does

1. **Cell QC and filtering**: load Cell Ranger matrices, compute QC metrics, flag doublets, optionally score ambient RNA, and save a filtered AnnData object.
2. **Guide processing and assignment**: parse guide reads, count UMIs, match guide references, assign single-guide / multi-guide / non-targeting cells, and summarize guide confidence.
3. **Batch diagnosis**: quantify donor, batch, condition, and guide confounding with PCA/UMAP summaries and recommended downstream handling.
4. **Escaped perturbation diagnosis**: support Mixscape-style KO / NP / NT classification and perturbation score summaries.
5. **Pseudobulk and differential expression**: aggregate by donor, condition, cell type, and guide; export edgeR-ready design matrices and contrasts.
6. **Pathway and visualization outputs**: volcano plots, heatmaps, enrichment tables, and a compact HTML report.

The defaults are intentionally conservative for large Perturb-seq/CROP-seq screens: preserve batch covariates for inference, aggregate to biological replicates before differential expression, and treat escaped perturbations as a first-class QC axis.

## Repository Layout

```text
config/                  Analysis configuration and sample/guide metadata templates
data/raw/                Cell Ranger inputs, FASTQs, and references
data/processed/          Stepwise pipeline outputs
docs/                    Architecture, reproducibility, methods, and glossary notes
notebooks/               Optional exploratory notebooks
scripts/                 R and utility scripts called by the workflow
src/perturbseq_pipeline/ Python package for reusable pipeline logic
tests/                   Unit tests for core data transformations
workflow/                Snakemake workflow and modular rules
results/final_reports/   Final summaries for collaborators or interview review
```

## Quick Start

```bash
conda env create -f environment.yml
conda activate perturbseq-pipeline
python -m perturbseq_pipeline make-demo --outdir data/demo
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --cores 4
```

For an HPC execution example:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --profile config/slurm
```

## Typical Inputs

- Cell Ranger filtered feature-barcode matrix directory or `.h5` matrix.
- Guide capture matrix or guide assignment table.
- `config/samples.csv`: sample, donor, batch, condition, and file paths.
- `config/guides_reference.csv`: guide IDs, target genes, non-targeting flags, and expected constructs.
- Optional `config/cell_annotations.yml`: curated cell-type labels or marker references.

## Key Outputs

- `data/processed/01_qc/adata_qc.h5ad`
- `data/processed/02_guides/adata_guides.h5ad`
- `data/processed/03_batch/batch_recommendations.md`
- `data/processed/03_mixscape/mixscape_scores.csv`
- `data/processed/04_pseudobulk/pseudobulk_counts.csv`
- `data/processed/05_de/de_results.csv`
- `data/processed/06_visualization/pathway_results.csv`
- `results/final_reports/summary_dashboard.html`

## Design Choices

- **Snakemake first** for inspectable, resumable, HPC-friendly workflow execution.
- **Python package modules** for testable data wrangling and QC logic.
- **R scripts** for edgeR/Mixscape-style methods that are common in genomics teams.
- **Config-driven analysis** so new screens can be run without editing code.
- **Audit-friendly outputs** at every step to make scientific decisions reviewable.

## Interview Talking Points

- How doublets, low-quality cells, and ambient RNA affect perturbation assignment.
- Why pseudobulk inference is preferred over naive cell-level tests for many replicated screens.
- How to detect and handle guide, donor, batch, and condition confounding.
- How escaped perturbations change interpretation of KO screens.
- How this repository could scale on SLURM or cloud batch execution.
- How processed perturbation matrices can feed downstream foundation-model training.

