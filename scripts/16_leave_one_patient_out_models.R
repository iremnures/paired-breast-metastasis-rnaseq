suppressPackageStartupMessages(library(DESeq2))

primary <- readRDS(
  "data/processed/deseq2_primary_dds.rds"
)

annotation <- read.csv(
  "data/processed/primary_liver_lung_counts_filtered.csv.gz",
  check.names = FALSE
)

previous_ranks <- readRDS(
  "data/processed/gsea_wald_rank_lists.rds"
)

raw <- counts(primary, normalized = FALSE)
metadata <- as.data.frame(colData(primary))
patients <- sort(unique(as.character(metadata$patient_code)))

stopifnot(
  length(patients) == 6,
  identical(colnames(raw), rownames(metadata)),
  !anyDuplicated(annotation$Geneid),
  all(table(metadata$patient_code, metadata$tissue) == 1L)
)

dir.create(
  "data/processed/leave_one_out",
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  "results/tables/leave_one_out",
  recursive = TRUE,
  showWarnings = FALSE
)

# Keep the primary ranks for the later common-universe comparison.
rank_lists <- list(primary = previous_ranks$primary)
summary_rows <- list()

for (patient in patients) {
  model_name <- paste0("no_", patient)

  cat("\n========================================\n")
  cat("Starting model:", model_name, "\n")
  cat("Time:", format(Sys.time()), "\n")

  keep <- as.character(metadata$patient_code) != patient
  meta <- metadata[keep, , drop = FALSE]

  # Remove the excluded patient's unused factor level.
  meta$patient_code <- droplevels(factor(meta$patient_code))
  meta$tissue <- factor(
    as.character(meta$tissue),
    levels = c("Liver", "Lung")
  )

  stopifnot(
    nrow(meta) == 10,
    nlevels(meta$patient_code) == 5,
    all(table(meta$patient_code, meta$tissue) == 1L)
  )

  # Reset inherited normalization before constructing each model.
  meta$sizeFactor <- NULL

  dds <- DESeqDataSetFromMatrix(
    countData = raw[, keep, drop = FALSE],
    colData = meta,
    design = ~ patient_code + tissue
  )

  stopifnot(
    is.null(sizeFactors(dds)),
    is.null(normalizationFactors(dds))
  )

  # Same defaults as the original scripts; sequential execution.
  dds <- DESeq(dds, parallel = FALSE)

  res <- results(
    dds,
    contrast = c("tissue", "Lung", "Liver"),
    alpha = 0.05
  )

  if (patient == "R18") {
    earlier <- readRDS("data/processed/deseq2_no_R18_dds.rds")
    earlier_res <- results(
      earlier,
      contrast = c("tissue", "Lung", "Liver"),
      alpha = 0.05
    )
    stopifnot(identical(rownames(res), rownames(earlier_res)))
    agreement <- all.equal(
      as.data.frame(res), as.data.frame(earlier_res),
      check.attributes = FALSE, tolerance = 1e-8
    )
    if (!isTRUE(agreement)) {
      print(agreement)
      stop("Corrected R18 model differs from the earlier R18 model.")
    }
    cat("\nCorrected R18 model matches the earlier R18 model.\n")
    rm(earlier, earlier_res)
  }

  tab <- data.frame(
    Geneid = rownames(res),
    baseMean = res$baseMean,
    log2FC_MLE = res$log2FoldChange,
    Wald_stat = res$stat,
    pvalue = res$pvalue,
    padj = res$padj
  )

  tab$gene_name <- annotation$gene_name[
    match(tab$Geneid, annotation$Geneid)
  ]
  tab$ensembl_id <- sub("\\.[0-9]+$", "", tab$Geneid)

  stopifnot(!anyDuplicated(tab$ensembl_id))

  write.csv(
    tab,
    paste0(
      "results/tables/leave_one_out/",
      model_name, "_gene_results.csv"
    ),
    row.names = FALSE
  )

  # Retain every finite Wald statistic, regardless of significance.
  finite <- is.finite(tab$Wald_stat)
  rank_tab <- tab[finite, , drop = FALSE]
  rank_tab <- rank_tab[
    order(-rank_tab$Wald_stat, rank_tab$ensembl_id),
    , drop = FALSE
  ]

  ranks <- setNames(
    rank_tab$Wald_stat,
    rank_tab$ensembl_id
  )
  rank_lists[[model_name]] <- ranks

  saveRDS(
    dds,
    paste0(
      "data/processed/leave_one_out/",
      model_name, "_dds.rds"
    )
  )

  saveRDS(
    ranks,
    paste0(
      "data/processed/leave_one_out/",
      model_name, "_wald_ranks.rds"
    )
  )

  fit_type <- attr(dispersionFunction(dds), "fitType")
  if (is.null(fit_type)) fit_type <- NA_character_

  summary_rows[[model_name]] <- data.frame(
    excluded_patient = patient,
    specimens = ncol(dds),
    genes_in_model = nrow(dds),
    finite_Wald_stats = sum(finite),
    nonfinite_Wald_stats = sum(!finite),
    genes_FDR_below_0.05 = sum(
      is.finite(tab$padj) & tab$padj < 0.05
    ),
    dispersion_fit = fit_type
  )

  # Update the summary after every completed model.
  write.csv(
    do.call(rbind, summary_rows),
    "results/tables/leave_one_out/model_summary.csv",
    row.names = FALSE
  )

  cat("\nCompleted:", model_name, "\n")
  print(summary_rows[[model_name]], row.names = FALSE)

  rm(dds, res)
  invisible(gc())
}

common_ids <- Reduce(
  intersect,
  lapply(rank_lists, names)
)

stopifnot(length(common_ids) > 0)

saveRDS(
  list(
    ranks = rank_lists,
    common_ids = common_ids,
    excluded_patients = patients
  ),
  "data/processed/leave_one_out/wald_rank_inputs.rds"
)

capture.output(
  sessionInfo(),
  file = "results/tables/leave_one_out/session_info.txt"
)

cat("\n--- All models completed ---\n")
print(do.call(rbind, summary_rows), row.names = FALSE)

cat("\n--- Genes available for GSEA ---\n")
print(lengths(rank_lists))
cat(
  "\nCommon finite genes across primary and all six models:",
  length(common_ids), "\n"
)
