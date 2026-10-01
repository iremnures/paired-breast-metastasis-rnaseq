# Paired transcriptomic comparison of liver and lung metastases in ER-positive breast cancer

## Overview

This analysis asks how transcriptomic profiles differ between
liver and lung metastasis specimens within the same patients.

The primary paired analysis identified 1,617 genes with
Benjamini–Hochberg adjusted p-values below 0.05: 1,212 with
higher estimated expression in Lung and 405 in Liver.

Hallmark enrichment identified 36 significant gene sets.
Of these, 24 retained the primary enrichment direction and
FDR < 0.05 in every leave-one-patient-out analysis.

The results support metastatic-site-associated differences
in bulk gene expression. They do not establish tumor-cell-
intrinsic adaptation independently of tissue composition.

## Dataset and analysis design

The source dataset is GSE316391:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316391

The primary cohort comprises six patients with matched liver
and lung metastasis specimens: R18, R23, R33, R36, R49 and R8.
There are 12 biological specimens.

RNA-seq library replicates were assessed for concordance and
raw counts were summed at the patient–tissue specimen level.
Library replicates were not treated as independent biological
replicates.

Genes were retained if they had counts of at least 10 in at
least six biological specimens. This rule was specified before
differential-expression results were inspected. The retained
analysis universe contained 25,667 genes.

DESeq2 version 1.50.2 was used with the patient-blocked design:

    ~ patient_code + tissue

The contrast was Lung versus Liver. Positive log2 fold changes
indicate Lung direction; negative values indicate Liver direction.

The clinical source labels liver specimens as "Liver/Bile Duct";
the analysis uses the shorter label "Liver".

## Sample-level quality control

VST-based PCA and sample distances showed both patient-specific
variation and site-associated structure.

R18LIV showed particularly strong separation along PC1.
It was flagged for inspection and retained in the primary
analysis because no independent technical QC evidence justified
its exclusion.

Samples were not removed because of PCA position or disagreement
with an expected biological result.

## Primary differential-expression results

All 25,667 retained genes had finite primary p-values and
available adjusted p-values.

At gene-level FDR < 0.05:

| Direction | Significant genes |
|---|---:|
| Lung | 1,212 |
| Liver | 405 |
| Total | 1,617 |

This imbalance in gene counts does not imply globally higher
transcription in Lung specimens or identify the cellular source
of the differences.

Paired expression inspection supported consistent organ-associated
patterns, including Liver-direction FGB and ALB expression and
Lung-direction SFTPB and SFTPC expression.

Some genes showed substantial patient heterogeneity, including
NKX2-1. UGT2A3 had positive raw counts in all liver specimens
and zero counts in all lung specimens; COL6A5 showed the
opposite pattern.

Observed zero counts do not establish biological absence.
Exact fold-change magnitudes for genes with all-zero counts
at one site require particular caution.

## Gene-level sensitivity to R18 exclusion

The model was refitted after removing both R18 specimens,
retaining the primary prefiltered gene universe and
re-estimating normalization and dispersion parameters.

The analysis excluding R18 identified 1,364 significant genes,
with 1,038 significant in both analyses.

All 1,617 primary-significant genes retained their estimated
effect direction. Within this selected gene set, the Spearman
correlation of unshrunken log2 fold changes was 0.963.
This is not a genome-wide correlation.

The seven organ-associated markers inspected in that analysis
retained direction and FDR significance.

## LFC shrinkage diagnostics

Comparison of default apeglm and random-start optimization
identified 12 genes with an absolute difference greater than
0.5 log2 units. This was a diagnostic threshold, not a
biological exclusion criterion.

The flagged genes were:

GALNT13, UGT1A5, CYP4F62P, PROK2, FGB, GRIA1, GABRB2,
CLDN2, SLN, IQSEC3-AS1, IGF1 and KCNJ16.

A third optimization method, nbinomR, did not consistently agree
with either method. UGT1A5 had unavailable posterior SD and
optimizer diagnostics under nbinomR.

These genes remain in the primary Wald-test results. Their
shrunken estimates are not used for definitive effect-size
prioritization while numerical disagreement remains unresolved.

Optimization convergence codes and agreement between methods
were not treated as proof of a globally optimal solution.
Estimates were not selectively combined across methods.

## Hallmark enrichment

Human MSigDB Hallmark version 2026.1.Hs was analyzed with
fgseaMultilevel using signed Wald statistics.

All 25,667 finite-statistic genes were ranked without filtering
on significance. The same gene universe and membership of
50 eligible gene sets were used across models.

BH adjustment was applied separately across the 50 gene sets
within each model. Pathway-level FDR is distinct from gene-level
FDR. Enrichment did not use the disputed shrinkage estimates.

The primary analysis identified 36 significant gene sets:
15 Lung-direction and 21 Liver-direction.

Selected prominent results were:

| Gene set | Primary NES | Direction |
|---|---:|---|
| Oxidative phosphorylation | -3.449 | Liver |
| MYC targets V1 | -3.196 | Liver |
| mTORC1 signaling | -2.563 | Liver |
| Inflammatory response | 2.315 | Lung |
| Epithelial–mesenchymal transition | 2.311 | Lung |
| Interferon gamma response | 2.056 | Lung |

NES is not an expression fold change or a direct measurement
of pathway activity.

## Leave-one-patient-out Hallmark sensitivity

Each patient's two specimens were removed together.
Six DESeq2 models were refitted using the original defaults
and the same paired design.

Every model retained 25,667 finite Wald statistics.
No enrichment NES or adjusted p-values were missing.

Among the 36 primary-significant gene sets:

| Sensitivity outcome | Gene sets |
|---|---:|
| Same direction and FDR < 0.05 in all six models | 24 |
| Same direction in all six, variable FDR support | 10 |
| Direction changed in at least one model | 2 |

Oxidative phosphorylation, MYC targets V1, inflammatory
response and EMT retained direction and FDR < 0.05 in all
six sensitivity models.

Coagulation changed from Liver to Lung direction after
excluding R18 or R8. Neither positive result met FDR < 0.05.
Its enrichment curves showed patterns at both ends of the
ranking; a pathway direction change did not imply reversal
of every member gene.

Heme metabolism changed direction after excluding R18.
It did not meet FDR < 0.05 after excluding R18, R33 or R8.

TGF-beta signaling retained Lung direction in all six models,
but met FDR < 0.05 only after excluding R8.

Complement did not meet FDR < 0.05 in the primary analysis.
It met the threshold only after excluding R18 or R8 and is
therefore treated as a sensitivity finding.

These models are overlapping subsets of the same dataset,
not independent validations. Removing a patient also reduces
the sample size. NES differences between models were not
formally tested.

Approximately 0.04% tied ranking statistics were reported
in the enrichment runs.

## Leading-edge inspection

For the primary versus R18-excluded comparison, leading-edge
Jaccard overlaps were approximately:

- Oxidative phosphorylation: 0.89
- MYC targets V1: 0.79
- Inflammatory response: 0.87
- EMT: 0.74

These describe gene-membership similarity, not significance
or biological activity.

Coagulation had zero leading-edge overlap when its selected
enrichment direction changed. Several strongly Liver-direction
genes, including FGA, nevertheless retained negative Wald
statistics.

Gene sets can share genes and should not be interpreted as
independent biological discoveries.

## Tissue-composition interpretation

Clinical records were matched using study, container and
internal case identifiers.

Tumor content and necrosis were unknown for 10 of 12 specimens.
Only R49 had records at both sites: 75%-99% tumor content
and zero recorded necrosis.

Tumor-content metadata were insufficient for covariate
adjustment. A shared interval does not establish equal tumor
fractions.

Exploratory marker inspection found:

- CD79A higher in Lung in 6/6 pairs; PTPRC, CD3D and LST1 in 5/6.
- DCN higher in Lung in 6/6 pairs; COL1A2 and LUM in 5/6,
  and COL1A1 in 4/6.
- PECAM1, VWF and CDH5 higher in Lung in 5/6 pairs.
- KRT18 higher in Liver in 6/6 pairs; EPCAM and KRT8 in 5/6.
- ALB higher in Liver in 6/6 pairs; APOA1 in 4/6.
- SFTPB and SFTPC higher in Lung in 6/6 pairs;
  SCGB1A1 in 5/6.

Immune, stromal/ECM and endothelial expression patterns are
compatible with composition and microenvironment contributions
to the site-associated signatures. They do not quantify cell
fractions or establish cellular origin.

Marker heatmap colours are within-gene z-scores. Descriptive
paired log-count differences are not DESeq2 model LFC estimates.

## Conclusion and remaining limitations

Matched liver and lung metastasis specimens showed substantial
site-associated transcriptomic differences. Major metabolic,
MYC-target, inflammatory and EMT-associated signatures persisted
across every single-patient exclusion tested.

The six-patient cohort, incomplete tumor-content metadata and
bulk measurements limit generalization and cellular attribution.
Independent validation and cell-resolved or spatial evidence
would be needed to establish tumor-cell-specific mechanisms.

The primary analysis retains all six patients. Optimizer-sensitive
shrinkage estimates remain an explicitly unresolved limitation.

## Supporting figures and records

- results/figures/selected_genes_paired_expression.pdf
- results/figures/ma_before_after_shrinkage.pdf
- results/figures/hallmark_enrichment_curves.pdf
- results/figures/hallmark_primary_vs_no_R18_summary.pdf
- results/figures/hallmark_leave_one_out_heatmap.pdf
- results/figures/exploratory_marker_heatmap.pdf
- results/tables/apeglm_results_with_diagnostic_flags.csv
- results/tables/leave_one_out/hallmark_stability_summary.csv
- results/tables/tumor_content_metadata_audit.csv
- docs/analysis_plan.md
- docs/qc_decisions.md
- docs/hallmark_gsea_findings.md
- docs/cell_composition_interpretation.md

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
