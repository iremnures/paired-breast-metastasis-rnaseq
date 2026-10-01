suppressPackageStartupMessages(library(DESeq2))

dds <- readRDS(
  "data/processed/deseq2_primary_dds.rds"
)

comparison <- read.csv(
  "results/tables/apeglm_three_method_comparison.csv",
  check.names = FALSE
)

# Diagnostic flags; these are not biological exclusion rules.
comparison$optimizer_sensitive <- (
  is.finite(comparison$optimizer_difference) &
  abs(comparison$optimizer_difference) > 0.5
)

comparison$R_posterior_SD_missing <-
  !is.finite(comparison$R_posterior_SD)

write.csv(
  comparison,
  "results/tables/apeglm_results_with_diagnostic_flags.csv",
  row.names = FALSE
)

# Genes chosen to inspect organ markers, patient influence,
# optimizer disagreement and unusually extreme estimates.
genes <- c(
  "APCS", "FGA", "FGB", "SFTPB",
  "SFTPC", "NKX2-1", "C8A", "SOSTDC1",
  "PROK2", "UGT1A5", "UGT2A3", "COL6A5"
)

idx <- match(genes, comparison$gene_name)
stopifnot(!anyNA(idx))

gene_ids <- comparison$Geneid[idx]
stopifnot(all(gene_ids %in% rownames(dds)))

norm <- counts(dds, normalized = TRUE)
meta <- as.data.frame(colData(dds))

patients <- unique(as.character(meta$patient_code))

# Verify that every patient has exactly one specimen per site.
pairs <- table(meta$patient_code, meta$tissue)
stopifnot(
  all(c("Liver", "Lung") %in% colnames(pairs)),
  all(pairs[, c("Liver", "Lung"), drop = FALSE] == 1L)
)

patient_colors <- setNames(
  hcl.colors(length(patients), palette = "Dark 3"),
  patients
)
patient_colors["R18"] <- "#B22222"

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save the plotted values for reproducibility.
expression_table <- do.call(
  rbind,
  lapply(seq_along(genes), function(i) {
    values <- as.numeric(norm[gene_ids[i], ])
    data.frame(
      Geneid = gene_ids[i],
      gene_name = genes[i],
      specimen_id = colnames(dds),
      patient_code = meta$patient_code,
      tissue = meta$tissue,
      normalized_count = values,
      log2_normalized_count_plus1 = log2(values + 1)
    )
  })
)

write.csv(
  expression_table,
  "results/tables/selected_genes_paired_expression.csv",
  row.names = FALSE
)

draw_plots <- function() {
  par(
    mfrow = c(3, 4),
    mar = c(3.2, 4.2, 2.5, 1),
    oma = c(3, 0, 3, 0)
  )

  for (i in seq_along(genes)) {
    values <- log2(as.numeric(norm[gene_ids[i], ]) + 1)
    upper <- max(values) + 0.8

    plot(
      NA,
      xlim = c(0.8, 2.2),
      ylim = c(0, upper),
      xaxt = "n",
      xlab = "",
      ylab = "log2(normalized count + 1)",
      main = genes[i],
      bty = "l"
    )
    axis(1, at = c(1, 2), labels = c("Liver", "Lung"))

    for (patient in patients) {
      liver <- which(
        as.character(meta$patient_code) == patient &
        as.character(meta$tissue) == "Liver"
      )
      lung <- which(
        as.character(meta$patient_code) == patient &
        as.character(meta$tissue) == "Lung"
      )

      lines(
        c(1, 2),
        values[c(liver, lung)],
        type = "b",
        pch = 16,
        col = patient_colors[patient],
        lwd = if (patient == "R18") 2.5 else 1.2
      )
    }
  }

  mtext(
    "Matched liver and lung specimens",
    outer = TRUE, side = 3, line = 1, font = 2
  )
  mtext(
    "Each line represents one patient; R18 is highlighted in red. Panel y-scales differ.",
    outer = TRUE, side = 1, line = 1
  )
}

pdf(
  "results/figures/selected_genes_paired_expression.pdf",
  width = 14,
  height = 10
)
draw_plots()
dev.off()

png(
  "results/figures/selected_genes_paired_expression.png",
  width = 2800,
  height = 2000,
  res = 200
)
draw_plots()
dev.off()

cat("\nPaired gene plots and diagnostic flags saved.\n")
