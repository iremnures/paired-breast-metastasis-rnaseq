suppressPackageStartupMessages(library(DESeq2))

annotation <- read.csv(
  "data/processed/primary_liver_lung_counts_filtered.csv.gz",
  check.names = FALSE
)

stopifnot(!anyDuplicated(annotation$Geneid))

model_files <- c(
  primary = "data/processed/deseq2_primary_dds.rds",
  no_R18 = "data/processed/deseq2_no_R18_dds.rds"
)

rank_lists <- list()

for (model_name in names(model_files)) {

  dds <- readRDS(model_files[[model_name]])

  res <- results(
    dds,
    contrast = c("tissue", "Lung", "Liver"),
    alpha = 0.05
  )

  tab <- data.frame(
    Geneid = rownames(res),
    baseMean = res$baseMean,
    log2FC_MLE = res$log2FoldChange,
    Wald_stat = res$stat,
    padj = res$padj
  )

  tab$gene_name <- annotation$gene_name[
    match(tab$Geneid, annotation$Geneid)
  ]

  # Remove Ensembl version suffixes, if present.
  tab$ensembl_id <- sub("\\.[0-9]+$", "", tab$Geneid)

  # Do not filter on significance.
  valid <- (
    is.finite(tab$Wald_stat) &
    !is.na(tab$ensembl_id) &
    nzchar(tab$ensembl_id)
  )

  cat("\n---", model_name, "---\n")
  cat("Genes before rank checks:", nrow(tab), "\n")
  cat("Excluded by rank checks:", sum(!valid), "\n")

  tab <- tab[valid, ]

  # Duplicate IDs require an explicit decision; do not choose
  # whichever statistic is most extreme.
  stopifnot(!anyDuplicated(tab$ensembl_id))

  tab <- tab[
    order(-tab$Wald_stat, tab$ensembl_id),
  ]

  ranks <- setNames(tab$Wald_stat, tab$ensembl_id)
  rank_lists[[model_name]] <- ranks

  write.csv(
    tab,
    paste0(
      "results/tables/gsea_wald_ranks_",
      model_name,
      ".csv"
    ),
    row.names = FALSE
  )

  cat("Genes retained:", length(ranks), "\n")
  cat("Positive Wald statistics:", sum(ranks > 0), "\n")
  cat("Negative Wald statistics:", sum(ranks < 0), "\n")
  cat("Zero Wald statistics:", sum(ranks == 0), "\n")
  cat("Repeated statistic values:", sum(duplicated(ranks)), "\n")

  cat("\nTop 5 Lung-direction genes:\n")
  print(
    head(tab[, c("gene_name", "Wald_stat", "padj")], 5),
    row.names = FALSE
  )

  cat("\nTop 5 Liver-direction genes:\n")
  print(
    tail(tab[, c("gene_name", "Wald_stat", "padj")], 5),
    row.names = FALSE
  )
}

saveRDS(
  rank_lists,
  "data/processed/gsea_wald_rank_lists.rds"
)

cat("\n--- Required packages ---\n")
cat(
  "fgsea installed:",
  requireNamespace("fgsea", quietly = TRUE),
  "\n"
)
cat(
  "msigdbr installed:",
  requireNamespace("msigdbr", quietly = TRUE),
  "\n"
)
