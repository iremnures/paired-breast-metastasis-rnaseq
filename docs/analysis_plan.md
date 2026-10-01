# M07 Bulk RNA-seq Analysis Plan

## Dataset
GSE316391

## Biological system
Metastatic ER-positive breast cancer.

## Primary question
How does the transcriptomic profile differ between liver and lung metastases within the same patients?

## Biological unit
Patient-tissue specimen.

RNA-seq library replicates from the same patient-tissue specimen are not treated as independent biological replicates.

## Replicate handling
RNA-seq library replicates were assessed for concordance and raw counts were summed at the biological specimen level.

## Primary cohort
Six patients with matched liver and lung metastases.

Total:
- 6 patients
- 6 liver specimens
- 6 lung specimens
- 12 biological specimens

## Primary model
Paired / patient-blocked design:

~ patient + tissue

## Primary contrast
Lung vs Liver

Interpretation:
- log2FC > 0: higher expression in Lung metastasis
- log2FC < 0: higher expression in Liver metastasis

## Low-count filtering
Genes are retained if they have:

count >= 10 in at least 6 biological specimens.

This filtering rule was specified before differential-expression results were inspected.

## Multiple testing
Benjamini-Hochberg false discovery rate correction.

## QC principles
Samples will not be removed because they fail to support an expected biological result.

Potential outliers will be evaluated using sample-level QC, PCA, sample distances, count distributions, and metadata.

## Interpretation boundary
Statistical significance alone will not be treated as evidence of biological importance.

Effect size, uncertainty, adjusted p-values, replicate-level/sample-level patterns, and study design will be considered together.

## Deviations
Any deviation from this analysis plan will be documented separately with justification.
