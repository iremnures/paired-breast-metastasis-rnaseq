tab <- read.csv(
  "results/tables/exploratory_marker_expression.csv",
  check.names = FALSE
)

genes <- unique(tab$gene_name)
paired_tables <- list()
summary_tables <- list()

for (gene in genes) {
  g <- tab[tab$gene_name == gene, ]
  liver <- g[g$tissue == "Liver", ]
  lung <- g[g$tissue == "Lung", ]

  stopifnot(
    nrow(liver) == 6,
    nrow(lung) == 6,
    !anyDuplicated(liver$patient_code),
    !anyDuplicated(lung$patient_code),
    setequal(liver$patient_code, lung$patient_code)
  )

  lung <- lung[match(liver$patient_code, lung$patient_code), ]

  difference <- (
    lung$log2_normalized_count_plus1 -
    liver$log2_normalized_count_plus1
  )

  paired_tables[[gene]] <- data.frame(
    gene_name = gene,
    group = liver$group,
    patient_code = liver$patient_code,
    Liver_raw_count = liver$raw_count,
    Lung_raw_count = lung$raw_count,
    Liver_normalized_count = liver$normalized_count,
    Lung_normalized_count = lung$normalized_count,
    Lung_minus_Liver_log2_count_plus1 = difference
  )

  summary_tables[[gene]] <- data.frame(
    gene_name = gene,
    group = g$group[1],
    Lung_higher_pairs = sum(difference > 0),
    Liver_higher_pairs = sum(difference < 0),
    equal_pairs = sum(difference == 0),
    median_paired_log_difference = median(difference),
    minimum_raw_count = min(g$raw_count),
    maximum_raw_count = max(g$raw_count),
    zero_count_specimens = sum(g$raw_count == 0)
  )
}

paired <- do.call(rbind, paired_tables)
summary <- do.call(rbind, summary_tables)
rownames(summary) <- NULL

write.csv(
  paired,
  "results/tables/exploratory_marker_paired_values.csv",
  row.names = FALSE
)

write.csv(
  summary,
  "results/tables/exploratory_marker_paired_summary.csv",
  row.names = FALSE
)

cat("\n--- Descriptive marker summary ---\n")
print(summary, row.names = FALSE)
