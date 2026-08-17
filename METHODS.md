# Methods and Decision Logic

The organizing idea is that every biological decision is represented as data, and
non-identifiable inference fails closed.

## Problem and Contract

Start in `README.md`, `config/analysis.yaml`, and one small input configuration. The bundled
`config/demo_h5ad.yaml` and `config/demo_csv.yaml` exercise matched representations of the same
screen. The input contract supports count CSVs, a combined Cell Ranger H5, a combined-feature
H5AD, a canonical H5AD, or a manifest of multiple runs. Manifest mode prefixes raw
10x barcodes with `sample_id`; every adapter emits the same canonical sparse AnnData object with
RNA counts in `X`/`layers['counts']`, guide counts in `obsm['guide_counts']`, and validated metadata
in `obs`. Missing sample ID, donor, batch, condition, or cell type fails before QC. Low-MOI
assumptions, guide thresholds, batch covariates, Mixscape settings, replicate definitions, and
references are configuration rather than hidden notebook state.

## Architecture

Open `workflow/Snakefile`. Snakemake is the orchestration and provenance layer;
`src/perturbseq_pipeline` contains unit-testable Python transformations; `workflow/scripts` contains
the Seurat, edgeR, and fgsea entry points. Stepwise outputs are inspectable CSV, JSON, PDF, and HTML
artifacts, and a failed workflow resumes from its last valid boundary.

Trace the DAG: input normalization → cell calling/QC and guide calling → RNA and guide-capture
batch diagnosis plus Mixscape → identifiability gate → replicate pseudobulk → edgeR → fgsea/report.

## Cell and Guide Decisions

Open `qc.py`, `guides.py`, and the `merge_metadata` rule. Filtered Cell Ranger/H5AD
inputs retain their upstream cell calls; custom demultiplexed matrices can use audited fixed
thresholds. QC adds library size, detected genes, mitochondrial fraction, and Scrublet.

Guide calling uses both minimum UMI support and top-guide fraction. With low MOI, two supported
guides are a multiplet even when one dominates. The merged `include_primary_analysis` flag requires
a called, QC-passing transcriptomic singlet with exactly one confident guide. Every excluded
barcode remains in the audit table with a reason.

Guide-capture batch diagnosis uses samples, not cells, as independent units. It tests assignment
rates, UMI and confidence summaries, and per-guide recovery across batches, applies BH correction,
and requires both statistical and effect-size thresholds. Missing replication is labeled rather
than converted into a cell-level p-value.

## Batch, Escape, and the Stop Rule

Open `expression_diagnostics.py`, `guide_diagnostics.py`, `design_gate.py`, and
`workflow/scripts/run_mixscape.R`. Expression batch diagnosis measures design sparsity and
association with uncorrected expression PCs; guide batch diagnosis separately tests replicate-level
capture and assignment metrics. Integration is visualization-only.
Mixscape assigns KO, NP (escaped/non-perturbed), and NT labels from local signatures.

The hard gate examines batch versus donor, assigned guide, target, and escape status. It represents each observed
contingency table as a bipartite graph. Multiple connected components mean exact aliasing, so the
effects cannot be separated. The JSON decision is an input to both pseudobulk branches;
`require_design_gate` raises before aggregation when blocked. This distinguishes a strong but
estimable batch effect from a completely confounded design.

## Replicate-Aware Inference

Open `pseudobulk.py` and `workflow/scripts/run_edger.R`. Counts are summed within sample ID, donor, batch,
condition, cell type, and target. Biological samples—not cells—are the units of inference. edgeR fits
quasi-likelihood target-versus-NT models within cell type, records under-replicated contrasts as
skipped, and refuses rank-deficient designs without silently dropping covariates.

The primary branch is intention-to-perturb and retains KO plus NP guide singlets. The sensitivity
branch keeps KO and NT only. A stronger KO-only effect suggests dilution by escaped perturbations;
an effect that is unstable across branches should be reported as such.

## Reproduction, Evidence, and Limitations

The reproducibility surface is `environment.yml`, `Makefile`,
`.github/workflows/ci.yml`, one run manifest, and `tests`. A validated run is:

```bash
make setup
make test
make dryrun-h5ad
make run
```

`make dryrun-csv` and `make run-csv` select the CSV adapter. Snakemake merges the selected input
configuration with `config/analysis.yaml` and writes the exact resolved configuration into the
processed-data directory before analysis. CI runs lint/unit tests and the full Conda-backed H5AD
workflow. Each run records source and configuration checksums, Git state, and core package
versions in `results/provenance/run_metadata.json`; exact R versions are pinned in
`environment.yml`. Provenance excludes timestamps, usernames, host details, locale, and local
installation paths. The deterministic fixture contains 12 sample-level
replicates, Poisson guide loading, and escaped perturbations, exercising the scientific branches
without protected data.

Close with limitations: production thresholds must be reviewed by sample; raw droplet matrices
need an ambient-aware upstream caller; Mixscape needs sufficient NT and per-target cells; exact
confounding cannot be repaired computationally; and gene sets/reference annotations must be
study-specific. These boundaries are visible rather than hidden.
