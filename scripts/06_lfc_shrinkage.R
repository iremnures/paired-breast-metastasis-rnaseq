suppressPackageStartupMessages(library(DESeq2))

stopifnot(requireNamespace("apeglm", quietly = TRUE))

# ------------------------------------------------------------
# 1. Load the primary fitted model
# ------------------------------------------------------------
dds <- readRDS(
  "data/processed/deseq2_primary_dds.rds"
)

coef_name <- "tissue_Lung_vs_Liver"

stopifnot(coef_name %in% resultsNames(dds))

# ------------------------------------------------------------
# 2. Original primary contrast
# Positive = higher in Lung; negative = higher in Liver
# ------------------------------------------------------------
res_raw <- results(
  dds,
  contrast = c("tissue", "Lung", "Liver"),
  alpha = 0.05
)

# ------------------------------------------------------------
# 3. Shrink the tissue coefficient
# Keep the original test results; no new effect-size threshold.
# ------------------------------------------------------------
res_shrunk <- lfcShrink(
  dds,
  coef = coef_name,
  res = res_raw,
  type = "apeglm",
  lfcThreshold = 0,
  svalue = FALSE,
  parallel = FALSE
)

stopifnot(
  identical(rownames(res_raw), rownames(res_shrunk))
)

# ------------------------------------------------------------
# 4. Annotation
# ------------------------------------------------------------
annotation <- read.csv(
  "data/processed/primary_liver_lung_counts_filtered.csv.gz",
  check.names = FALSE
)

stopifnot(!anyDuplicated(annotation$Geneid))

idx <- match(rownames(res_raw), annotation$Geneid)
stopifnot(!anyNA(idx))

# ------------------------------------------------------------
# 5. Side-by-side result table
# ------------------------------------------------------------
tab <- data.frame(
  Geneid = rownames(res_raw),
  gene_name = annotation$gene_name[idx],
  baseMean = res_raw$baseMean,
  log2FC_MLE = res_raw$log2FoldChange,
  lfcSE_MLE = res_raw$lfcSE,
  log2FC_shrunken = res_shrunk$log2FoldChange,
  posterior_SD = res_shrunk$lfcSE,
  pvalue = res_raw$pvalue,
  padj = res_raw$padj
)

tab$LFC_change <- tab$log2FC_shrunken - tab$log2FC_MLE

tab$abs_LFC_reduction <-
  abs(tab$log2FC_MLE) - abs(tab$log2FC_shrunken)

tab <- tab[order(tab$padj, na.last = TRUE), ]

sig <- tab[
  !is.na(tab$padj) & tab$padj < 0.05,
]

# ------------------------------------------------------------
# 6. Save results
# ------------------------------------------------------------
dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  res_shrunk,
  "data/processed/deseq2_primary_lfc_apeglm.rds"
)

write.csv(
  tab,
  "results/tables/deseq2_lung_vs_liver_apeglm_all_genes.csv",
  row.names = FALSE
)

write.csv(
  sig,
  "results/tables/deseq2_lung_vs_liver_apeglm_fdr05.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Console summary
# ------------------------------------------------------------
cat("\n--- LFC shrinkage summary ---\n")
cat("Genes in table:", nrow(tab), "\n")
cat("Primary Wald FDR < 0.05:", nrow(sig), "\n")

cat("\nTop 15 genes by primary adjusted p-value:\n")
print(
  head(
    sig[, c(
      "gene_name", "baseMean",
      "log2FC_MLE", "lfcSE_MLE",
      "log2FC_shrunken", "posterior_SD", "padj"
    )],
    15
  ),
  row.names = FALSE
)

# Organ markers and genes affected by R18 exclusion
inspect_genes <- c(
  "APCS", "CRP", "FGA", "GC",
  "SFTPA1", "SFTPB", "SFTPC",
  "C8A", "NKX2-1", "SOSTDC1", "SFTA3"
)

cat("\nSelected genes:\n")
selected <- tab[
  match(inspect_genes, tab$gene_name, nomatch = 0L),
  c(
    "gene_name", "baseMean",
    "log2FC_MLE", "lfcSE_MLE",
    "log2FC_shrunken", "posterior_SD", "padj"
  )
]
print(selected, row.names = FALSE)

finite <- (
  is.finite(sig$log2FC_MLE) &
  is.finite(sig$log2FC_shrunken)
)

largest_reductions <- sig[finite, ]
largest_reductions <- largest_reductions[
  order(
    largest_reductions$abs_LFC_reduction,
    decreasing = TRUE
  ),
]

cat("\nLargest absolute-LFC reductions among primary-significant genes:\n")
print(
  head(
    largest_reductions[, c(
      "gene_name", "baseMean",
      "log2FC_MLE", "lfcSE_MLE",
      "log2FC_shrunken", "posterior_SD",
      "abs_LFC_reduction", "padj"
    )],
    15
  ),
  row.names = FALSE
)
