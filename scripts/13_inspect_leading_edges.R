analysis <- readRDS(
  "data/processed/hallmark_gsea_analysis.rds"
)

primary_ranks <- read.csv(
  "results/tables/gsea_wald_ranks_primary.csv",
  check.names = FALSE
)

sensitivity_ranks <- read.csv(
  "results/tables/gsea_wald_ranks_no_R18.csv",
  check.names = FALSE
)

selected <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_COAGULATION",
  "HALLMARK_COMPLEMENT"
)

primary <- analysis$results$primary
sensitivity <- analysis$results$no_R18

stopifnot(
  all(selected %in% primary$pathway),
  all(selected %in% sensitivity$pathway),
  !anyDuplicated(primary_ranks$ensembl_id),
  !anyDuplicated(sensitivity_ranks$ensembl_id)
)

gene_tables <- list()
summaries <- list()

for (pathway in selected) {
  p <- match(pathway, primary$pathway)
  s <- match(pathway, sensitivity$pathway)

  le_primary <- unique(primary$leadingEdge[[p]])
  le_sensitivity <- unique(sensitivity$leadingEdge[[s]])
  ids <- union(le_primary, le_sensitivity)

  pi <- match(ids, primary_ranks$ensembl_id)
  si <- match(ids, sensitivity_ranks$ensembl_id)

  stopifnot(!anyNA(pi), !anyNA(si))

  tab <- data.frame(
    pathway = rep(pathway, length(ids)),
    ensembl_id = ids,
    gene_name = primary_ranks$gene_name[pi],
    leading_edge_primary = ids %in% le_primary,
    leading_edge_no_R18 = ids %in% le_sensitivity,
    primary_Wald_stat = primary_ranks$Wald_stat[pi],
    no_R18_Wald_stat = sensitivity_ranks$Wald_stat[si]
  )

  tab$leading_edge_both <- (
    tab$leading_edge_primary & tab$leading_edge_no_R18
  )

  # Display order only; this is not a GSEA contribution score.
  tab <- tab[
    order(-abs(tab$primary_Wald_stat), tab$ensembl_id),
  ]

  gene_tables[[pathway]] <- tab

  summaries[[pathway]] <- data.frame(
    pathway = pathway,
    primary_LE_genes = length(le_primary),
    no_R18_LE_genes = length(le_sensitivity),
    shared_LE_genes = length(intersect(le_primary, le_sensitivity)),
    union_LE_genes = length(ids),
    Jaccard_overlap = if (length(ids) > 0) {
      length(intersect(le_primary, le_sensitivity)) / length(ids)
    } else {
      NA_real_
    }
  )

  cat("\n---", pathway, "---\n")
  cat("First 10 genes by absolute primary Wald statistic:\n")
  print(
    head(tab[, c(
      "gene_name",
      "leading_edge_primary",
      "leading_edge_no_R18",
      "primary_Wald_stat",
      "no_R18_Wald_stat"
    )], 10),
    row.names = FALSE
  )
}

summary_table <- do.call(rbind, summaries)
rownames(summary_table) <- NULL

write.csv(
  do.call(rbind, gene_tables),
  "results/tables/selected_hallmark_leading_edge_genes.csv",
  row.names = FALSE
)

write.csv(
  summary_table,
  "results/tables/selected_hallmark_leading_edge_overlap.csv",
  row.names = FALSE
)

cat("\n--- Leading-edge overlap ---\n")
print(summary_table, row.names = FALSE)
