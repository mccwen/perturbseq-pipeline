# Methods Notes

## Cell QC

Cells are filtered using total counts, detected genes, mitochondrial fraction, and doublet score. In production, thresholds should be reviewed by sample and by donor rather than globally accepted without inspection.

## Guide Calling

Guide assignment uses UMI support and dominance of the top guide within each cell. Cells are labeled as `single`, `multi`, `ambiguous`, or `unassigned`. Non-targeting controls are preserved as their own class for perturbation modeling.

## Batch Diagnosis

The pipeline reports covariate balance across perturbations and recommends whether donor, batch, and condition should be included in the differential expression design. Batch integration is treated as visualization-only unless a downstream method explicitly supports corrected counts for inference.

## Escaped Perturbation Diagnosis

Mixscape-style labels separate non-targeting cells, cells with strong perturbation signatures, and cells that carry a guide but lack the expected expression shift. The primary analysis can include KO and NP singlets, while sensitivity analysis can restrict to KO-only cells.

## Differential Expression

The recommended production mode is pseudobulk aggregation followed by edgeR quasi-likelihood testing with donor and batch covariates. This reduces pseudoreplication risk and makes results easier to audit across biological replicates.

## Foundation Model Readiness

Stepwise outputs preserve perturbation labels, confidence metrics, covariates, cell states, escaped perturbation flags, and replicate-aware expression summaries. These are the metadata fields needed to convert a screen into supervised examples for perturbation foundation models.

