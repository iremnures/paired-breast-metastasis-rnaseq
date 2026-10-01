library(DESeq2)

# 1. Primary fitted object
primary <- readRDS(
  "data/processed/deseq2_primary_dds.rds"
)

stopifnot(
  all(c("R18LIV", "R18LUN") %in% colnames(primary)),
  all(c("patient_code", "tissue") %in% names(colData(primary)))
)

# 2. Remove the complete R18 pair
keep <- !colnames(primary) %in% c("R18LIV", "R18LUN")

meta <- as.data.frame(colData(primary))[keep, , drop = FALSE]
meta$patient_code <- droplevels(factor(meta$patient_code))
meta$tissue <- factor(
  as.character(meta$tissue),
  levels = c("Liver", "Lung")
)

# Keep metadata needed by the new model.
# Rebuild from raw counts to avoid reusing fitted estimates.
meta <- meta[, c("patient_code", "tissue"), drop = FALSE]

pair_table <- table(meta$patient_code, meta$tissue)

stopifnot(
  nrow(meta) == 10L,
  nrow(pair_table) == 5L,
  all(pair_table == 1L)
)

cat("\nRemaining patient_code pairs:\n")
print(pair_table)

sensitivity <- DESeqDataSetFromMatrix(
  countData = counts(primary, normalized = FALSE)[
    , keep, drop = FALSE
  ],
  colData = meta,
  design = ~ patient_code + tissue
)

# Same prefiltered gene universe; refit all estimates.
sensitivity <- DESeq(sensitivity)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  sensitivity,
  "data/processed/deseq2_no_R18_dds.rds"
)

# 3. Same contrast and FDR target in both models
contrast <- c("tissue", "Lung", "Liver")

res_primary <- results(
  primary,
  contrast = contrast,
  alpha = 0.05
)

res_sensitivity <- results(
  sensitivity,
  contrast = contrast,
  alpha = 0.05
)

stopifnot(
  identical(rownames(res_primary), rownames(res_sensitivity))
)

comparison <- data.frame(
  Geneid = rownames(res_primary),
  primary_log2FC = res_primary$log2FoldChange,
  primary_lfcSE = res_primary$lfcSE,
  primary_padj = res_primary$padj,
  no_R18_log2FC = res_sensitivity$log2FoldChange,
  no_R18_lfcSE = res_sensitivity$lfcSE,
  no_R18_padj = res_sensitivity$padj
)

# 4. Add annotation without changing row order
annotation <- read.csv(
  "data/processed/primary_liver_lung_counts_filtered.csv.gz",
  check.names = FALSE
)

stopifnot(!anyDuplicated(annotation$Geneid))

comparison$gene_name <- annotation$gene_name[
  match(comparison$Geneid, annotation$Geneid)
]

comparison$delta_log2FC <-
  comparison$no_R18_log2FC - comparison$primary_log2FC

comparison$same_direction <-
  sign(comparison$primary_log2FC) ==
  sign(comparison$no_R18_log2FC)

primary_sig <- (
  !is.na(comparison$primary_padj) &
  comparison$primary_padj < 0.05
)

sensitivity_sig <- (
  !is.na(comparison$no_R18_padj) &
  comparison$no_R18_padj < 0.05
)

finite_lfc <- (
  is.finite(comparison$primary_log2FC) &
  is.finite(comparison$no_R18_log2FC)
)

# 5. Save complete comparison
comparison <- comparison[, c(
  "Geneid", "gene_name",
  "primary_log2FC", "primary_lfcSE", "primary_padj",
  "no_R18_log2FC", "no_R18_lfcSE", "no_R18_padj",
  "delta_log2FC", "same_direction"
)]

write.csv(
  comparison,
  "results/tables/deseq2_primary_vs_no_R18.csv",
  row.names = FALSE
)

# 6. Summary
cat("\n--- R18 sensitivity summary ---\n")
cat("Primary FDR < 0.05:", sum(primary_sig), "\n")
cat("Without R18 FDR < 0.05:", sum(sensitivity_sig), "\n")
cat(
  "Significant in both models:",
  sum(primary_sig & sensitivity_sig), "\n"
)
cat(
  "Primary-significant genes with an estimable sensitivity LFC:",
  sum(primary_sig & finite_lfc), "\n"
)
cat(
  "Same direction among these genes:",
  sum(comparison$same_direction[primary_sig & finite_lfc]),
  "\n"
)
cat(
  "Spearman LFC correlation among primary-significant genes:",
  round(
    cor(
      comparison$primary_log2FC[primary_sig & finite_lfc],
      comparison$no_R18_log2FC[primary_sig & finite_lfc],
      method = "spearman"
    ),
    3
  ),
  "\n"
)

# 7. Organ-marker comparison
markers <- c("APCS", "CRP", "FGA", "GC", "SFTPA1", "SFTPB", "SFTPC")

cat("\nOrgan-marker comparison:\n")
print(
  comparison[
    comparison$gene_name %in% markers,
    c(
      "gene_name",
      "primary_log2FC", "no_R18_log2FC",
      "primary_padj", "no_R18_padj"
    )
  ],
  row.names = FALSE
)

# 8. Largest changes among primary-significant genes
changed <- comparison[primary_sig & finite_lfc, ]
changed <- changed[
  order(abs(changed$delta_log2FC), decreasing = TRUE),
]

cat("\nLargest LFC changes among primary-significant genes:\n")
print(
  head(changed, 15),
  row.names = FALSE
)
