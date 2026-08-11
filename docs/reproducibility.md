# Reproducibility

## Environment

Use `environment.yml` for a full conda environment, or `pyproject.toml` for Python-only development.

```bash
conda env create -f environment.yml
conda activate perturbseq-pipeline
```

## Workflow

All production steps should be launched through Snakemake:

```bash
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --cores 8
```

## Versioning

- Keep raw data immutable under `data/raw`.
- Store exact references and guide annotations under `data/references` or `config`.
- Commit configs, workflow code, docs, and tests.
- Do not commit large `.h5ad`, FASTQ, or Cell Ranger matrix outputs.

## Review Checklist

- QC thresholds justified and sample-aware.
- Guide assignment summary reviewed before DE.
- Donor, batch, condition, and guide confounding checked.
- Pseudobulk design matrix matches the biological question.
- Escaped perturbations reported and handled consistently.
- Final report includes limitations and sensitivity analyses.

