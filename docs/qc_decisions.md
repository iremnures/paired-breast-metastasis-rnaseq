# QC Decisions and Exploratory Observations

## 07.1 Sample-level QC

VST-transformed sample-to-sample distances and PCA indicated structured variation among the paired liver and lung metastasis specimens.

In the paired PCA, liver and lung metastases from several patients showed a broadly consistent directional separation, suggesting a potential metastatic-site-associated transcriptomic signal.

At the same time, substantial patient-specific variation was present. R18LIV showed particularly strong separation along PC1 and was therefore flagged for further inspection rather than automatically excluded.

No sample was removed based on PCA position or because it failed to support an expected biological result.

These observations are exploratory only. Formal testing of the liver-versus-lung difference will be performed using a patient-blocked DESeq2 model.

## Initial differential-expression interpretation check

The primary paired DESeq2 analysis identified strong metastatic-site-associated expression differences.

Several canonical organ-associated genes showed consistent paired directionality across all six patients. Liver-associated markers including APCS, CRP, FGA and GC were consistently higher in liver metastasis specimens, whereas lung-associated surfactant genes including SFTPA1, SFTPB and SFTPC were consistently higher in lung metastasis specimens.

These patterns indicate that the bulk RNA-seq contrast captures a strong metastatic-site signal. However, this signal cannot be interpreted as exclusively tumor-cell-intrinsic because bulk tissue composition and the local organ microenvironment may contribute to the observed expression differences.

Tumor percentage metadata were unavailable for most specimens and therefore cannot be included as an adjustment variable.

R18LIV remains in the primary analysis because no independent technical QC evidence currently justifies exclusion. A sensitivity analysis excluding the R18 matched pair will be performed later.

Primary interpretation:
site-associated transcriptomic differences between matched liver and lung metastases.

Tumor-cell-intrinsic adaptation will not be claimed from bulk RNA-seq alone.

## Sensitivity analysis excluding the R18 matched pair

The paired DESeq2 model was refitted after excluding both R18LIV
and R18LUN, retaining the primary prefiltered gene universe and
re-estimating normalization and dispersion parameters.

At FDR < 0.05, the primary analysis identified 1,617 genes and the
sensitivity analysis identified 1,364; 1,038 were significant in both.

All 1,617 primary-significant genes retained their estimated effect
direction. Spearman correlation of unshrunken log2 fold changes
within this selected gene set was 0.963.

All seven inspected organ-associated markers retained their
direction and FDR significance after R18 exclusion.

These findings support a site-associated signal that persists
without R18. Some gene-level effect estimates and significance
statuses remain sensitive to inclusion of this patient.

R18 remains in the primary analysis. Bulk tissue composition and
microenvironment contributions remain unresolved.

## LFC shrinkage diagnostics and zero-count patterns

Comparison of default apeglm optimization and random-start
optimization identified 12 genes with an absolute difference
greater than 0.5 log2 units. This threshold is diagnostic only.

A third optimization method, nbinomR, did not uniformly agree
with either method. UGT1A5 had unavailable posterior SD and
optimizer diagnostics under nbinomR.

These genes remain in the primary Wald-test results but are
flagged as optimizer-sensitive. Their shrunken effect estimates
will not be used for definitive effect-size prioritization until
the numerical disagreement is resolved.

Raw-count inspection showed:
- UGT2A3: positive counts in all six liver specimens and zero
  counts in all six lung specimens.
- COL6A5: zero counts in all six liver specimens and positive
  counts in all six lung specimens.
- NKX2-1: heterogeneous expression across patient pairs,
  including zero counts at both sites for R36.
- FGB: higher normalized expression in liver for all six pairs.
- UGT1A5: positive raw counts in all specimens.

Zero observed counts do not establish biological absence.
Genes with all-zero counts at one site retain their primary
test results, but exact fold-change magnitudes require caution.

Bulk tissue composition remains an unresolved interpretation
limitation.

## Normalization correction and verification

During validation, script 16 was found to inherit primary-model
size factors through specimen metadata. The script was corrected
to remove inherited size factors before constructing each
patient-exclusion model.

Each model now re-estimates normalization from its retained
specimens. Assertions check that normalization factors are
initially absent and that the corrected R18-excluded model
matches the earlier independently refitted R18 model.

Stages 16-18 and their outputs were regenerated. The corrected
results retain the main sensitivity summary: 24 primary-significant
Hallmark sets retain direction and FDR < 0.05 in all six models;
10 retain direction with variable FDR support; two change direction.

The corrected stages were also rerun in the separate reproduction
directory and compared with the corrected project results.
Earlier uncorrected LOO outputs are superseded.
