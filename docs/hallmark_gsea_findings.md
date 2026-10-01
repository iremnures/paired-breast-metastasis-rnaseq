# Hallmark GSEA findings

## Analysis approach

Human MSigDB Hallmark gene sets (version 2026.1.Hs) were tested
using fgseaMultilevel and signed DESeq2 Wald statistics for
Lung versus Liver.

Both models used the same 25,667 genes and identical membership
for all 50 eligible gene sets. No significance filtering was
applied to the ranked gene lists.

Positive NES indicates Lung-direction enrichment.
Negative NES indicates Liver-direction enrichment.

BH correction was applied separately across the 50 gene sets
within each model. GSEA did not use the optimizer-sensitive
apeglm effect estimates.

## Primary and R18 sensitivity results

The primary model identified 36 gene sets with FDR < 0.05:
15 Lung-direction and 21 Liver-direction.

The model excluding R18 identified 31 gene sets with FDR < 0.05:
15 Lung-direction and 16 Liver-direction.

Thirty gene sets met FDR < 0.05 in both models, with consistent
NES direction. Six met the threshold only in the primary model,
one only in the model excluding R18, and 13 in neither model.

## Selected findings retained after excluding R18

- Oxidative phosphorylation: Liver direction; NES -3.449 to -2.879.
- MYC targets V1: Liver direction; NES -3.196 to -3.212.
- Inflammatory response: Lung direction; NES 2.315 to 2.640.
- Epithelial-mesenchymal transition: Lung direction;
  NES 2.311 to 2.454.

All four gene sets had FDR < 0.05 in both models.
Their leading-edge Jaccard overlaps were approximately
0.89, 0.79, 0.87 and 0.74, respectively.

Enrichment curves were consistent with these directions.

## Findings sensitive to excluding R18

Adipogenesis, coagulation, cholesterol homeostasis, reactive
oxygen species pathway, TGF-beta signaling and heme metabolism
met FDR < 0.05 only in the primary model.

Coagulation and heme metabolism changed NES direction.
Neither had FDR < 0.05 after excluding R18.

Coagulation retained strongly negative Wald statistics for
several genes, including FGA. Its enrichment curves showed
patterns at both ends of the ranking, with a change in the
dominant signed excursion. Leading-edge overlap was zero;
this does not imply that all member genes reversed direction.

Complement retained Lung-direction enrichment, with NES
increasing from 1.307 to 1.858 and FDR changing from 0.0522
to approximately 2.34e-6. Leading-edge Jaccard overlap was 0.73.

## Interpretation limits

These are site-associated expression signatures in bulk tissue.
They do not establish tumor-cell-specific pathway activation,
causality, or independence from tissue composition.

Excluding R18 also reduces the number of patient pairs from
six to five. Changes in FDR cannot be attributed solely to
R18 biology or interpreted automatically as loss of a signal.

NES is not an expression fold change. Differences between
model NES values were not formally tested.

Gene sets can share genes and are not independent biological
findings. Leading-edge overlap is descriptive, not a
significance test.

Tied ranking statistics affected approximately 0.04% of each
ranked list. This warning was retained as a diagnostic note.

Sensitivity to removing the other patients has not yet been
assessed.

## Outputs

- results/tables/hallmark_sensitivity_summary.csv
- results/tables/selected_hallmark_leading_edge_genes.csv
- results/tables/selected_hallmark_leading_edge_overlap.csv
- results/figures/hallmark_enrichment_curves.pdf
- results/figures/hallmark_primary_vs_no_R18_summary.pdf
- results/tables/hallmark_gsea_session_info.txt

## Leave-one-patient-out sensitivity analysis

Each patient's liver and lung specimens were removed together.
Six models were refitted from counts using the original
DESeq2 defaults and the design ~ patient_code + tissue.

All models retained 25,667 finite Wald statistics. GSEA used
the same ranked gene universe and the same 50 Hallmark gene
sets, with BH correction applied separately within each model.

Of the 36 primary-significant gene sets:
- 24 retained the primary NES direction and FDR < 0.05 in all
  six leave-one-patient-out models.
- 10 retained direction in all six models, with variable
  FDR-threshold support.
- 2 changed direction in at least one model: coagulation and
  heme metabolism.

Oxidative phosphorylation, MYC targets V1, inflammatory
response and epithelial-mesenchymal transition retained
direction and FDR < 0.05 in all six sensitivity models.

Coagulation changed from Liver to Lung direction after
excluding R18 or R8; neither positive result met FDR < 0.05.

Heme metabolism changed direction after excluding R18 and
did not meet FDR < 0.05 after excluding R18, R33 or R8.

E2F targets and G2M checkpoint retained Liver direction but
did not meet FDR < 0.05 after excluding R23.

TNFA signaling via NFKB retained Lung direction but did not
meet FDR < 0.05 after excluding R36.

TGF-beta signaling retained Lung direction in all six models,
but met FDR < 0.05 only after excluding R8.

Complement, which did not meet FDR < 0.05 in the primary model,
met the threshold only after excluding R18 or R8.

These overlapping subsets assess sensitivity within the
available dataset; they are not independent validations.
No patient was excluded from the primary analysis on the
basis of these results.

No NES or adjusted p-values were missing. Approximately 0.04%
tied ranking statistics were reported in each sensitivity run.

Additional outputs:
- results/tables/leave_one_out/hallmark_stability_summary.csv
- results/tables/leave_one_out/hallmark_sensitivity_changes.csv
- results/figures/hallmark_leave_one_out_heatmap.pdf

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
