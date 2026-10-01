suppressPackageStartupMessages(library(DESeq2))

dds <- readRDS(
  "data/processed/deseq2_primary_dds.rds"
)

original <- readRDS(
  "data/processed/deseq2_primary_lfc_apeglm.rds"
)

raw <- results(
  dds,
  contrast = c("tissue", "Lung", "Liver"),
  alpha = 0.05
)

# Keep the entire gene universe so prior adaptation uses
# the same data as in the original shrinkage analysis.
set.seed(20260930)

alternative <- lfcShrink(
  dds,
  coef = "tissue_Lung_vs_Liver",
  res = raw,
  type = "apeglm",
  apeMethod = "nbinomC*",
  parallel = FALSE
)

stopifnot(
  identical(rownames(raw), rownames(original)),
  identical(rownames(raw), rownames(alternative))
)

annotation <- read.csv(
  "data/processed/primary_liver_lung_counts_filtered.csv.gz",
  check.names = FALSE
)

stopifnot(!anyDuplicated(annotation$Geneid))

idx <- match(rownames(raw), annotation$Geneid)
stopifnot(!anyNA(idx))

tab <- data.frame(
  Geneid = rownames(raw),
  gene_name = annotation$gene_name[idx],
  log2FC_MLE = raw$log2FoldChange,
  default_apeglm = original$log2FoldChange,
  random_start_apeglm = alternative$log2FoldChange,
  primary_padj = raw$padj
)

tab$optimizer_difference <-
  tab$random_start_apeglm - tab$default_apeglm

write.csv(
  tab,
  "results/tables/apeglm_optimizer_comparison.csv",
  row.names = FALSE
)

saveRDS(
  alternative,
  "data/processed/deseq2_primary_apeglm_random_start_check.rds"
)

cat("\n--- Package versions ---\n")
cat("DESeq2:", as.character(packageVersion("DESeq2")), "\n")
cat("apeglm:", as.character(packageVersion("apeglm")), "\n")

cat("\n--- Original prior information ---\n")
print(priorInfo(original))

genes <- c(
  "APCS", "SFTPB", "SFTPC",
  "FGB", "NKX2-1", "SOSTDC1", "SFTA3", "C8A",
  "UGT2A3", "COL6A5"
)

cat("\n--- Selected genes: optimizer comparison ---\n")
print(
  tab[match(genes, tab$gene_name, nomatch = 0L), ],
  row.names = FALSE
)

finite <- is.finite(tab$optimizer_difference)
sig <- !is.na(tab$primary_padj) & tab$primary_padj < 0.05

# Diagnostic threshold only, not a biological selection rule.
cat(
  "\nGenes with optimizer difference > 0.5 log2 units:",
  sum(abs(tab$optimizer_difference[finite]) > 0.5),
  "\n"
)

cat(
  "Among primary-significant genes:",
  sum(abs(tab$optimizer_difference[finite & sig]) > 0.5),
  "\n"
)
