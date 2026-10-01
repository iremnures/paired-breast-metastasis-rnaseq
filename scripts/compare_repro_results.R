suppressPackageStartupMessages(library(DESeq2))

run_dir <- trimws(readLines(
  "repro_checks/latest_run.txt",
  warn = FALSE
)[1])

stopifnot(
  dir.exists(run_dir),
  file.exists(file.path(run_dir, "execution_completed.txt"))
)

checks <- list()

add_check <- function(name, passed) {
  checks[[name]] <<- data.frame(
    check = name,
    passed = isTRUE(passed)
  )
}

cat("\n--- Reproduction directory ---\n")
cat(run_dir, "\n")

input_files <- c(
  "counts_by_biological_specimen.csv.gz",
  "specimen_metadata.csv",
  "primary_liver_lung_counts.csv.gz",
  "primary_liver_lung_metadata.csv",
  "primary_liver_lung_counts_filtered.csv.gz"
)

for (filename in input_files) {
  original <- read.csv(
    file.path("data/processed", filename),
    check.names = FALSE
  )
  reproduced <- read.csv(
    file.path(run_dir, "data/processed", filename),
    check.names = FALSE
  )

  add_check(
    paste("Input:", filename),
    isTRUE(all.equal(
      original, reproduced,
      check.attributes = FALSE,
      tolerance = 0
    ))
  )
}

model_files <- c(
  primary = "data/processed/deseq2_primary_dds.rds",
  no_R18 = "data/processed/leave_one_out/no_R18_dds.rds",
  no_R23 = "data/processed/leave_one_out/no_R23_dds.rds",
  no_R33 = "data/processed/leave_one_out/no_R33_dds.rds",
  no_R36 = "data/processed/leave_one_out/no_R36_dds.rds",
  no_R49 = "data/processed/leave_one_out/no_R49_dds.rds",
  no_R8 = "data/processed/leave_one_out/no_R8_dds.rds"
)

model_summaries <- list()

for (model in names(model_files)) {
  original <- readRDS(model_files[[model]])
  reproduced <- readRDS(file.path(
    run_dir, model_files[[model]]
  ))

  stopifnot(
    identical(rownames(original), rownames(reproduced)),
    identical(colnames(original), colnames(reproduced))
  )

  a <- as.data.frame(results(
    original,
    contrast = c("tissue", "Lung", "Liver"),
    alpha = 0.05
  ))
  b <- as.data.frame(results(
    reproduced,
    contrast = c("tissue", "Lung", "Liver"),
    alpha = 0.05
  ))

  add_check(
    paste("DESeq2 results:", model),
    isTRUE(all.equal(
      a, b,
      check.attributes = FALSE,
      tolerance = 1e-8
    ))
  )

  add_check(
    paste("Model counts:", model),
    isTRUE(all.equal(
      counts(original), counts(reproduced),
      tolerance = 0
    ))
  )

  model_summaries[[model]] <- data.frame(
    model = model,
    original_significant_genes = sum(
      is.finite(a$padj) & a$padj < 0.05
    ),
    reproduced_significant_genes = sum(
      is.finite(b$padj) & b$padj < 0.05
    )
  )

  rm(original, reproduced)
  invisible(gc())
}

old_gsea <- readRDS(
  "data/processed/leave_one_out/hallmark_analysis.rds"
)
new_gsea <- readRDS(file.path(
  run_dir,
  "data/processed/leave_one_out/hallmark_analysis.rds"
))

gsea_details <- list()

for (model in names(old_gsea$results)) {
  a <- old_gsea$results[[model]]
  b <- new_gsea$results[[model]]

  stopifnot(
    setequal(a$pathway, b$pathway),
    !anyDuplicated(a$pathway),
    !anyDuplicated(b$pathway)
  )

  b <- b[match(a$pathway, b$pathway), ]

  gsea_details[[model]] <- data.frame(
    model = model,
    pathway = a$pathway,
    original_NES = a$NES,
    reproduced_NES = b$NES,
    absolute_NES_difference = abs(a$NES - b$NES),
    original_padj = a$padj,
    reproduced_padj = b$padj,
    same_direction = sign(a$NES) == sign(b$NES),
    same_FDR_status = (
      (is.finite(a$padj) & a$padj < 0.05) ==
      (is.finite(b$padj) & b$padj < 0.05)
    )
  )
}

gsea <- do.call(rbind, gsea_details)
rownames(gsea) <- NULL

gsea_summary <- do.call(rbind, lapply(
  split(gsea, gsea$model),
  function(tab) {
    data.frame(
      model = tab$model[1],
      pathways = nrow(tab),
      direction_changes = sum(
        !tab$same_direction, na.rm = TRUE
      ),
      FDR_status_changes = sum(!tab$same_FDR_status),
      missing_NES_comparisons = sum(
        !is.finite(tab$absolute_NES_difference)
      ),
      maximum_NES_difference = max(
        tab$absolute_NES_difference,
        na.rm = TRUE
      )
    )
  }
))

old_flags <- read.csv(
  "results/tables/apeglm_results_with_diagnostic_flags.csv",
  check.names = FALSE
)
new_flags <- read.csv(file.path(
  run_dir,
  "results/tables/apeglm_results_with_diagnostic_flags.csv"
), check.names = FALSE)

old_sensitive <- old_flags$Geneid[
  old_flags$optimizer_sensitive %in% TRUE
]
new_sensitive <- new_flags$Geneid[
  new_flags$optimizer_sensitive %in% TRUE
]

add_check(
  "Same optimizer-sensitive gene set",
  setequal(old_sensitive, new_sensitive)
)

check_table <- do.call(rbind, checks)
rownames(check_table) <- NULL
model_summary <- do.call(rbind, model_summaries)

write.csv(
  check_table,
  "docs/reproducibility_comparison.csv",
  row.names = FALSE
)
write.csv(
  model_summary,
  "docs/reproduced_gene_result_counts.csv",
  row.names = FALSE
)
write.csv(
  gsea,
  "results/tables/reproduced_gsea_comparison.csv",
  row.names = FALSE
)

cat("\n--- Input and model comparisons ---\n")
print(check_table, row.names = FALSE)

cat("\n--- Significant gene counts ---\n")
print(model_summary, row.names = FALSE)

cat("\n--- GSEA comparison ---\n")
print(gsea_summary, row.names = FALSE)

cat("\n--- GSEA direction or FDR-status changes ---\n")
changed <- (
  is.na(gsea$same_direction) |
  !gsea$same_direction |
  !gsea$same_FDR_status
)
if (any(changed)) {
  print(gsea[changed, ], row.names = FALSE)
} else {
  cat("None.\n")
}

cat("\n--- Optimizer-sensitive genes ---\n")
cat("Original:", length(old_sensitive), "\n")
cat("Reproduced:", length(new_sensitive), "\n")
cat("Only in original:\n")
print(setdiff(old_sensitive, new_sensitive))
cat("Only in reproduced:\n")
print(setdiff(new_sensitive, old_sensitive))
