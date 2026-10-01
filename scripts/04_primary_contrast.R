library(DESeq2)

# ============================================================
# M07 - Primary DESeq2 contrast
# Lung vs Liver metastases
# ============================================================

# ------------------------------------------------------------
# 1. Load fitted DESeq2 object
# ------------------------------------------------------------

dds <- readRDS(
  "data/processed/deseq2_primary_dds.rds"
)

# ------------------------------------------------------------
# 2. Extract primary contrast
# Positive LFC = higher in Lung
# Negative LFC = higher in Liver
# ------------------------------------------------------------

res <- results(
  dds,
  contrast = c(
    "tissue",
    "Lung",
    "Liver"
  ),
  alpha = 0.05
)

# ------------------------------------------------------------
# 3. Convert to data frame
# ------------------------------------------------------------

res_df <- as.data.frame(res)

res_df$Geneid <- rownames(res_df)

# ------------------------------------------------------------
# 4. Add gene annotation
# ------------------------------------------------------------

annotation <- read.csv(
  "data/processed/primary_liver_lung_counts_filtered.csv.gz",
  check.names = FALSE
)

annotation <- annotation[
  ,
  c(
    "Geneid",
    "gene_name",
    "Chr",
    "Start",
    "End",
    "Strand",
    "Length"
  )
]

res_df <- merge(
  res_df,
  annotation,
  by = "Geneid",
  all.x = TRUE,
  sort = FALSE
)

# ------------------------------------------------------------
# 5. Reorder useful columns
# ------------------------------------------------------------

res_df <- res_df[
  ,
  c(
    "Geneid",
    "gene_name",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "Chr",
    "Start",
    "End",
    "Strand",
    "Length"
  )
]

# ------------------------------------------------------------
# 6. Sort by adjusted p-value
# ------------------------------------------------------------

res_ordered <- res_df[
  order(
    res_df$padj,
    na.last = TRUE
  ),
]

# ------------------------------------------------------------
# 7. Save complete result table
# ------------------------------------------------------------

write.csv(
  res_ordered,
  "results/tables/deseq2_lung_vs_liver_all_genes.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. FDR-significant subset
# ------------------------------------------------------------

sig <- res_ordered[
  !is.na(res_ordered$padj)
  & res_ordered$padj < 0.05,
]

write.csv(
  sig,
  "results/tables/deseq2_lung_vs_liver_fdr05.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Direction counts
# ------------------------------------------------------------

sig_lung <- sig[
  sig$log2FoldChange > 0,
]

sig_liver <- sig[
  sig$log2FoldChange < 0,
]

# ------------------------------------------------------------
# 10. Print summary
# ------------------------------------------------------------

cat("\nDESeq2 results summary:\n")

print(
  summary(res)
)

cat(
  "\nGenes tested:",
  nrow(res_df),
  "\n"
)

cat(
  "Genes with padj < 0.05:",
  nrow(sig),
  "\n"
)

cat(
  "FDR-significant, higher in Lung:",
  nrow(sig_lung),
  "\n"
)

cat(
  "FDR-significant, higher in Liver:",
  nrow(sig_liver),
  "\n"
)

cat(
  "\nGenes with NA adjusted p-value:",
  sum(is.na(res_df$padj)),
  "\n"
)

cat(
  "\nTop 20 genes by adjusted p-value:\n"
)

print(
  res_ordered[
    1:min(
      20,
      nrow(res_ordered)
    ),
    c(
      "Geneid",
      "gene_name",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "pvalue",
      "padj"
    )
  ]
)
