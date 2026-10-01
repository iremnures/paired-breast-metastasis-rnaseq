suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
})

dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load ranks and use a common finite gene universe
# ------------------------------------------------------------
rank_lists <- readRDS(
  "data/processed/gsea_wald_rank_lists.rds"
)

stopifnot(all(c("primary", "no_R18") %in% names(rank_lists)))

common_ids <- intersect(
  names(rank_lists$primary),
  names(rank_lists$no_R18)
)

stopifnot(length(common_ids) > 0)

cat("\n--- Rank universe ---\n")
cat("Primary original:", length(rank_lists$primary), "\n")
cat("No R18 original:", length(rank_lists$no_R18), "\n")
cat("Common genes:", length(common_ids), "\n")

for (model in c("primary", "no_R18")) {
  ranks <- rank_lists[[model]][common_ids]

  stopifnot(
    all(is.finite(ranks)),
    !anyDuplicated(names(ranks))
  )

  ranks <- ranks[order(-ranks, names(ranks))]
  rank_lists[[model]] <- ranks
}

# ------------------------------------------------------------
# 2. Retrieve and cache human Hallmark gene sets
# ------------------------------------------------------------
cache <- "data/metadata/msigdb_hallmark_human.rds"

if (file.exists(cache)) {
  hallmark <- readRDS(cache)
} else {
  hallmark <- msigdbr(
    species = "Homo sapiens",
    collection = "H"
  )
  saveRDS(hallmark, cache)
}

stopifnot(
  all(c("gs_name", "ensembl_gene") %in% names(hallmark))
)

valid <- (
  !is.na(hallmark$ensembl_gene) &
  nzchar(hallmark$ensembl_gene)
)

mapping <- unique(data.frame(
  pathway = hallmark$gs_name[valid],
  ensembl_id = sub(
    "\\.[0-9]+$", "",
    hallmark$ensembl_gene[valid]
  )
))

write.csv(
  mapping,
  "data/metadata/hallmark_ensembl_mapping.csv",
  row.names = FALSE
)

pathways <- lapply(
  split(mapping$ensembl_id, mapping$pathway),
  unique
)

# Identical pathway membership in both models.
pathways_common <- lapply(
  pathways,
  function(ids) intersect(ids, common_ids)
)

coverage <- data.frame(
  pathway = names(pathways),
  mapped_genes = lengths(pathways),
  genes_in_rank_universe = lengths(pathways_common)
)

coverage$eligible <- (
  coverage$genes_in_rank_universe >= 15 &
  coverage$genes_in_rank_universe <= 500
)

write.csv(
  coverage,
  "results/tables/hallmark_gene_coverage.csv",
  row.names = FALSE
)

pathways_used <- pathways_common[coverage$eligible]
stopifnot(length(pathways_used) > 0)

cat("\n--- Hallmark coverage ---\n")
cat("Available gene sets:", length(pathways), "\n")
cat("Eligible gene sets:", length(pathways_used), "\n")

if ("db_version" %in% names(hallmark)) {
  cat("MSigDB version:", unique(hallmark$db_version), "\n")
}

# ------------------------------------------------------------
# 3. Run both models
# BH correction is applied separately within each model.
# ------------------------------------------------------------
results_list <- list()

for (model in c("primary", "no_R18")) {

  set.seed(20260930)

  result <- fgseaMultilevel(
    pathways = pathways_used,
    stats = rank_lists[[model]],
    minSize = 15,
    maxSize = 500,
    eps = 0,
    scoreType = "std",
    BPPARAM = BiocParallel::SerialParam()
  )

  result <- as.data.frame(result)
  result <- result[order(result$padj, na.last = TRUE), ]
  results_list[[model]] <- result

  # Preserve list columns in RDS; flatten for CSV.
  csv_result <- result
  csv_result$leadingEdge <- vapply(
    csv_result$leadingEdge,
    paste,
    collapse = ";",
    FUN.VALUE = character(1)
  )

  write.csv(
    csv_result,
    paste0("results/tables/hallmark_gsea_", model, ".csv"),
    row.names = FALSE
  )

  sig <- !is.na(result$padj) & result$padj < 0.05

  cat("\n---", model, "Hallmark results ---\n")
  cat("Gene sets tested:", nrow(result), "\n")
  cat("Missing adjusted p-values:", sum(is.na(result$padj)), "\n")
  cat("FDR < 0.05:", sum(sig), "\n")
  cat("Lung-direction:", sum(sig & result$NES > 0), "\n")
  cat("Liver-direction:", sum(sig & result$NES < 0), "\n")

  cat("\nTop 15 by adjusted p-value:\n")
  print(
    head(
      result[, c("pathway", "NES", "pval", "padj", "size")],
      15
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 4. Compare primary and sensitivity results
# ------------------------------------------------------------
primary <- results_list$primary[, c("pathway", "NES", "padj")]
names(primary)[2:3] <- c("primary_NES", "primary_padj")

sensitivity <- results_list$no_R18[, c("pathway", "NES", "padj")]
names(sensitivity)[2:3] <- c("no_R18_NES", "no_R18_padj")

comparison <- merge(
  primary,
  sensitivity,
  by = "pathway",
  all = TRUE
)

comparison$same_direction <- (
  sign(comparison$primary_NES) ==
  sign(comparison$no_R18_NES)
)

comparison <- comparison[
  order(comparison$primary_padj, na.last = TRUE),
]

write.csv(
  comparison,
  "results/tables/hallmark_primary_vs_no_R18.csv",
  row.names = FALSE
)

saveRDS(
  list(
    ranks = rank_lists,
    pathways = pathways_used,
    results = results_list,
    comparison = comparison,
    msigdb_version = if ("db_version" %in% names(hallmark)) {
      unique(hallmark$db_version)
    } else {
      NA_character_
    }
  ),
  "data/processed/hallmark_gsea_analysis.rds"
)

capture.output(
  sessionInfo(),
  file = "results/tables/hallmark_gsea_session_info.txt"
)

cat("\n--- Primary vs no R18: top 15 primary gene sets ---\n")
print(head(comparison, 15), row.names = FALSE)
