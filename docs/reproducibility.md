# Reproducibility verification

## Scope

The numbered analysis workflow was executed in a separate
working directory using source counts, the clinical workbook,
and the cached MSigDB Hallmark 2026.1.Hs definitions.

After correcting inherited normalization factors in script 16,
stages 16-18 were rerun in that directory. Unchanged upstream
outputs came from the preceding full workflow execution.

Verification used the same installed R and package environment.
It does not establish portability to a different environment.

## Verified comparisons

- All five reconstructed input tables matched exactly.
- Counts matched for the primary and six patient-exclusion models.
- DESeq2 result tables matched within all.equal tolerance 1e-8.
- Significant gene counts matched in all seven models.
- All 350 Hallmark NES values matched exactly.
- Hallmark direction and FDR-threshold status matched throughout.
- The same 12 optimizer-sensitive genes were identified.

GSEA p-values and posterior shrinkage estimates were not required
to match exactly by this comparison.

## Significant gene counts at FDR < 0.05

| Model | Genes |
|---|---:|
| Primary | 1,617 |
| Excluding R18 | 1,364 |
| Excluding R23 | 818 |
| Excluding R33 | 1,469 |
| Excluding R36 | 1,018 |
| Excluding R49 | 1,721 |
| Excluding R8 | 1,624 |

## Evidence

- docs/reproducibility_comparison.csv
- docs/reproduced_gene_result_counts.csv
- results/tables/reproduced_gsea_comparison.csv
- docs/software_versions.csv
- docs/environment_session_info.txt

Execution logs are retained locally under repro_checks/.

## Remaining limitation

Reproduction of optimizer-sensitive flags does not resolve the
underlying shrinkage disagreement. Those estimates remain
unsuitable for definitive effect-size prioritization.
