.PHONY: help setup lint test demo dryrun dryrun-csv dryrun-h5ad run run-csv run-h5ad clean

INPUT_CONFIG ?= config/demo_h5ad.yaml

# Keep globally installed Python packages out of the project environment.
export PYTHONNOUSERSITE := 1
export XDG_CACHE_HOME := $(CURDIR)/.cache
export NUMBA_CACHE_DIR := $(CURDIR)/.cache/numba
export MPLCONFIGDIR := $(CURDIR)/.cache/matplotlib

help:
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "%-12s %s\n", $$1, $$2}'

setup: ## Reproduce the validated Conda and R environment
	conda env update -f environment.yml --prune
	conda run -n perturbseq-pipeline bash -c 'PATH="$$CONDA_PREFIX/bin:$$PATH" Rscript workflow/scripts/install_r_packages.R'

lint: ## Run Python lint checks
	conda run -n perturbseq-pipeline ruff check src tests

test: ## Run unit tests
	conda run -n perturbseq-pipeline python -m pytest -q

demo: ## Regenerate the synthetic low-MOI dataset
	conda run -n perturbseq-pipeline python -m perturbseq_pipeline make-demo --outdir demo

dryrun: ## Validate the DAG for INPUT_CONFIG (H5AD by default)
	conda run -n perturbseq-pipeline bash -c 'PATH="$$CONDA_PREFIX/bin:$$PATH" snakemake --snakefile workflow/Snakefile --configfile $(INPUT_CONFIG) --cores 1 --dry-run'

dryrun-csv: ## Validate the CSV input workflow
	$(MAKE) dryrun INPUT_CONFIG=config/demo_csv.yaml

dryrun-h5ad: ## Validate the H5AD input workflow
	$(MAKE) dryrun INPUT_CONFIG=config/demo_h5ad.yaml

run: ## Run the workflow for INPUT_CONFIG (H5AD by default)
	conda run -n perturbseq-pipeline bash -c 'PATH="$$CONDA_PREFIX/bin:$$PATH" snakemake --snakefile workflow/Snakefile --configfile $(INPUT_CONFIG) --cores 4'

run-csv: ## Run the complete CSV input workflow
	$(MAKE) run INPUT_CONFIG=config/demo_csv.yaml

run-h5ad: ## Run the complete H5AD input workflow
	$(MAKE) run INPUT_CONFIG=config/demo_h5ad.yaml

clean: ## Remove processed outputs and reports
	rm -rf data/processed results .snakemake .cache .ruff_cache .pytest_cache \
		src/perturbseq_pipeline/__pycache__ tests/__pycache__
