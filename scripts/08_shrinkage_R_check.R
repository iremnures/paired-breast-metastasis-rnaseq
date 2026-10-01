suppressPackageStartupMessages(library(DESeq2))

dds <- readRDS(
  "data/processed/deseq2_primary_dds.rds"
)

comparison <- read.csv(
  "results/tables/apeglm_optimizer_comparison.csv",
  check.names = FALSE
)

raw <- results(
  dds,
  contrast = c("tissue", "Lung", "Liver"),
  alpha = 0.05
)

# Same model, gene universe and prior-adaptation input.
# Only the optimization method changes.
fit_R <- lfcShrink(
  dds,
  coef = "tissue_Lung_vs_Liver",
  res = raw,
  type = "apeglm",
  apeMethod = "nbinomR",
  returnList = TRUE,
  parallel = FALSE
)

res_R <- fit_R[[1]]
details_R <- fit_R[[2]]

idx <- match(comparison$Geneid, rownames(res_R))
stopifnot(!anyNA(idx))

comparison$R_apeglm <- res_R$log2FoldChange[idx]
comparison$R_posterior_SD <- res_R$lfcSE[idx]

comparison$R_minus_default <-
  comparison$R_apeglm - comparison$default_apeglm

comparison$R_minus_random <-
  comparison$R_apeglm - comparison$random_start_apeglm

saveRDS(
  fit_R,
  "data/processed/deseq2_primary_apeglm_R_check.rds"
)

write.csv(
  comparison,
  "results/tables/apeglm_three_method_comparison.csv",
  row.names = FALSE
)

flag <- (
  is.finite(comparison$optimizer_difference) &
  abs(comparison$optimizer_difference) > 0.5
)

flagged <- comparison[flag, ]

cat("\n--- Prior used by nbinomR ---\n")
print(priorInfo(res_R))

cat("\n--- Previously flagged genes: three methods ---\n")
print(
  flagged[, c(
    "gene_name", "log2FC_MLE",
    "default_apeglm", "random_start_apeglm",
    "R_apeglm", "R_posterior_SD",
    "primary_padj"
  )],
  row.names = FALSE
)

cat("\n--- R optimizer diagnostics for flagged genes ---\n")
if (!is.null(details_R$diag)) {
  diagnostic_idx <- match(flagged$Geneid, rownames(res_R))
  print(details_R$diag[diagnostic_idx, , drop = FALSE])
} else {
  cat("No diagnostic matrix returned.\n")
}

cat("\n--- FGB sample-level expression ---\n")
fgb_id <- comparison$Geneid[comparison$gene_name == "FGB"]
stopifnot(length(fgb_id) == 1L)

sample_table <- data.frame(
  specimen_id = colnames(dds),
  patient_code = colData(dds)$patient_code,
  tissue = colData(dds)$tissue,
  raw_count = as.numeric(counts(dds)[fgb_id, ]),
  normalized_count = as.numeric(
    counts(dds, normalized = TRUE)[fgb_id, ]
  )
)

print(sample_table, row.names = FALSE)
