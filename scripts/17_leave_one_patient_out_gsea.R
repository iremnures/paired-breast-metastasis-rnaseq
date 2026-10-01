suppressPackageStartupMessages(library(fgsea))

inputs <- readRDS(
  "data/processed/leave_one_out/wald_rank_inputs.rds"
)

previous <- readRDS(
  "data/processed/hallmark_gsea_analysis.rds"
)

common_ids <- inputs$common_ids
pathways <- previous$pathways

# Reuse the primary result only if its inputs are unchanged.
stopifnot(
  setequal(common_ids, names(previous$ranks$primary)),
  identical(
    inputs$ranks$primary[common_ids],
    previous$ranks$primary[common_ids]
  ),
  length(pathways) == 50,
  all(unlist(pathways, use.names = FALSE) %in% common_ids)
)

results_list <- list(primary = previous$results$primary)

flatten_result <- function(result) {
  tab <- as.data.frame(result)
  tab$leadingEdge <- vapply(
    tab$leadingEdge,
    paste,
    collapse = ";",
    FUN.VALUE = character(1)
  )
  tab
}

for (model in setdiff(names(inputs$ranks), "primary")) {
  cat("\n--- Starting GSEA:", model, "---\n")

  ranks <- inputs$ranks[[model]][common_ids]
  stopifnot(
    all(is.finite(ranks)),
    !anyDuplicated(names(ranks))
  )
  ranks <- ranks[order(-ranks, names(ranks))]

  set.seed(20260930)

  result <- fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = 15,
    maxSize = 500,
    eps = 0,
    scoreType = "std",
    BPPARAM = BiocParallel::SerialParam()
  )

  result <- as.data.frame(result)
  result <- result[order(result$padj, na.last = TRUE), ]
  results_list[[model]] <- result

  saveRDS(
    result,
    paste0(
      "data/processed/leave_one_out/",
      model, "_hallmark_gsea.rds"
    )
  )

  write.csv(
    flatten_result(result),
    paste0(
      "results/tables/leave_one_out/",
      model, "_hallmark_gsea.csv"
    ),
    row.names = FALSE
  )

  sig <- is.finite(result$padj) & result$padj < 0.05

  cat("Gene sets tested:", nrow(result), "\n")
  cat("Missing FDR values:", sum(!is.finite(result$padj)), "\n")
  cat("FDR < 0.05:", sum(sig), "\n")
  cat("Lung-direction:", sum(sig & result$NES > 0, na.rm = TRUE), "\n")
  cat("Liver-direction:", sum(sig & result$NES < 0, na.rm = TRUE), "\n")
}

primary <- results_list$primary
loo_models <- setdiff(names(results_list), "primary")

NES <- sapply(loo_models, function(model) {
  result <- results_list[[model]]
  idx <- match(primary$pathway, result$pathway)
  stopifnot(!anyNA(idx))
  result$NES[idx]
})

FDR <- sapply(loo_models, function(model) {
  result <- results_list[[model]]
  idx <- match(primary$pathway, result$pathway)
  result$padj[idx]
})

same_direction <- sweep(
  sign(NES),
  1,
  sign(primary$NES),
  FUN = "=="
)

summary_table <- data.frame(
  pathway = primary$pathway,
  primary_NES = primary$NES,
  primary_padj = primary$padj,
  primary_FDR_below_0.05 = (
    is.finite(primary$padj) & primary$padj < 0.05
  ),
  LOO_models_same_direction = rowSums(
    same_direction, na.rm = TRUE
  ),
  LOO_models_FDR_below_0.05 = rowSums(
    is.finite(FDR) & FDR < 0.05
  ),
  LOO_models_missing_NES = rowSums(!is.finite(NES)),
  LOO_models_missing_FDR = rowSums(!is.finite(FDR))
)

# Add each model's values so threshold changes remain inspectable.
for (model in loo_models) {
  summary_table[[paste0(model, "_NES")]] <- NES[, model]
  summary_table[[paste0(model, "_padj")]] <- FDR[, model]
}

summary_table <- summary_table[
  order(summary_table$primary_padj, na.last = TRUE),
]

write.csv(
  summary_table,
  "results/tables/leave_one_out/hallmark_stability_summary.csv",
  row.names = FALSE
)

saveRDS(
  list(
    results = results_list,
    summary = summary_table,
    pathways = pathways,
    common_ids = common_ids,
    msigdb_version = previous$msigdb_version
  ),
  "data/processed/leave_one_out/hallmark_analysis.rds"
)

capture.output(
  sessionInfo(),
  file = "results/tables/leave_one_out/gsea_session_info.txt"
)

cat("\n--- Primary-significant pathways: sensitivity summary ---\n")

columns <- c(
  "pathway",
  "primary_NES",
  "primary_padj",
  "LOO_models_same_direction",
  "LOO_models_FDR_below_0.05",
  "LOO_models_missing_NES",
  "LOO_models_missing_FDR"
)

print(
  summary_table[
    summary_table$primary_FDR_below_0.05,
    columns
  ],
  row.names = FALSE
)

cat("\n--- All pathways: number of LOO models with FDR < 0.05 ---\n")
print(table(summary_table$LOO_models_FDR_below_0.05))

cat("\n--- All pathways: number of LOO models retaining primary direction ---\n")
print(table(summary_table$LOO_models_same_direction))
