# Exploratory cell-composition interpretation

## Purpose

Cell- and organ-associated genes were inspected to contextualize
Hallmark enrichment in bulk metastatic specimens. This was a
descriptive analysis, not cell-fraction estimation or a new
significance test.

## Clinical metadata audit

All 12 specimens were matched to clinical records using study,
container and internal case identifiers.

Tumor content and necrosis were unknown for 10 of 12 specimens.
Only R49 had records at both sites: 75%-99% tumor content and
zero recorded necrosis. These ranges do not establish equal
tumor fractions. Tumor content was not added as a model covariate.

All matched specimens were labelled invasive ductal carcinoma.
The original clinical site label "Liver/Bile Duct" was retained
in the audit; the analysis uses the label "Liver".

## Descriptive expression patterns

Based on DESeq2-normalized counts within patient pairs:

- CD79A was higher in Lung in 6/6 pairs; PTPRC, CD3D and LST1
  were higher in Lung in 5/6 pairs.
- DCN was higher in Lung in 6/6 pairs; COL1A2 and LUM in 5/6,
  and COL1A1 in 4/6.
- PECAM1, VWF and CDH5 were higher in Lung in 5/6 pairs.
- KRT18 was higher in Liver in 6/6 pairs; EPCAM and KRT8 in 5/6.
- ALB was higher in Liver in 6/6 pairs; APOA1 in 4/6.
- SFTPB and SFTPC were higher in Lung in 6/6 pairs;
  SCGB1A1 in 5/6.

Immune, stromal and endothelial panel genes had positive raw
counts in all specimens. Some immune-marker counts were low.
SFTPC had one zero-count specimen and SCGB1A1 had two.

TTR was omitted from the panel because no matching gene was
found in the retained analysis universe. This does not establish
biological absence.

## Interpretation limits

Immune and stromal/ECM-associated expression patterns accompany
Lung-direction inflammatory and EMT enrichment. These findings
are compatible with tissue-composition and microenvironment
contributions but do not establish their magnitude or cellular
source.

Epithelial-marker expression does not establish tumor purity.
Organ-associated expression does not by itself establish
contamination or justify excluding a specimen.

Heatmap colours are within-gene z-scores. Paired differences
are differences in log2(normalized count + 1), not model LFCs.

## Outputs

- results/tables/tumor_content_metadata_audit.csv
- results/tables/exploratory_marker_paired_summary.csv
- results/tables/exploratory_marker_paired_values.csv
- results/figures/exploratory_marker_heatmap.pdf
