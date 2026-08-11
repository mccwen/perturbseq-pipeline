.PHONY: help lint test demo dryrun run clean

help:
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "%-12s %s\n", $$1, $$2}'

lint: ## Run Python lint checks
	ruff check src tests

test: ## Run unit tests
	pytest -q

demo: ## Create a tiny demo dataset
	python -m perturbseq_pipeline make-demo --outdir data/demo

dryrun: ## Validate Snakemake DAG
	snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --cores 1 --dry-run

run: ## Run the full configured workflow
	snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --cores 4

clean: ## Remove generated demo and processed outputs
	rm -rf data/demo data/processed/*/* results/final_reports/*

