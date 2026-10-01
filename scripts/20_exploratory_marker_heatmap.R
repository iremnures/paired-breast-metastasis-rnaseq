suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

dds <- readRDS("data/processed/deseq2_primary_dds.rds")

annotation <- read.csv(
  "data/processed/primary_liver_lung_counts_filtered.csv.gz",
  check.names = FALSE
)

panels <- list(
  "Immune-associated" = c("PTPRC", "CD3D", "CD79A", "LST1"),
  "Stroma / ECM-associated" = c("COL1A1", "COL1A2", "DCN", "LUM"),
  "Endothelial-associated" = c("PECAM1", "VWF", "CDH5"),
  "Epithelial-associated" = c("EPCAM", "KRT8", "KRT18"),
  "Liver-associated" = c("ALB", "APOA1"),
  "Lung-associated" = c("SFTPB", "SFTPC", "SCGB1A1")
)

genes <- unlist(panels, use.names = FALSE)
groups <- rep(names(panels), lengths(panels))

matches <- lapply(genes, function(gene) {
  which(
    !is.na(annotation$gene_name) &
    annotation$gene_name == gene &
    annotation$Geneid %in% rownames(dds)
  )
})

mapping_check <- data.frame(
  gene_name = genes,
  group = groups,
  matched_gene_ids = lengths(matches)
)

write.csv(
  mapping_check,
  "results/tables/exploratory_marker_mapping_check.csv",
  row.names = FALSE
)

if (any(lengths(matches) != 1L)) {
  print(mapping_check[lengths(matches) != 1L, ], row.names = FALSE)
  stop("Missing or ambiguous gene mapping; inspect before plotting.")
}

ids <- annotation$Geneid[unlist(matches)]
meta <- as.data.frame(colData(dds))

sample_order <- order(
  as.character(meta$patient_code),
  match(as.character(meta$tissue), c("Liver", "Lung"))
)

raw <- counts(dds)[ids, sample_order, drop = FALSE]
norm <- counts(dds, normalized = TRUE)[
  ids, sample_order, drop = FALSE
]
meta <- meta[sample_order, , drop = FALSE]

log_values <- log2(norm + 1)
gene_sd <- apply(log_values, 1, sd)

# Constant rows have no relative variation to display.
z <- sweep(log_values, 1, rowMeans(log_values), "-")
z <- sweep(z, 1, ifelse(gene_sd > 0, gene_sd, 1), "/")

long <- do.call(rbind, lapply(seq_along(genes), function(i) {
  data.frame(
    Geneid = ids[i],
    gene_name = genes[i],
    group = groups[i],
    specimen_id = colnames(norm),
    patient_code = meta$patient_code,
    tissue = meta$tissue,
    raw_count = as.numeric(raw[i, ]),
    normalized_count = as.numeric(norm[i, ]),
    log2_normalized_count_plus1 = as.numeric(log_values[i, ]),
    row_z_score = as.numeric(z[i, ]),
    constant_expression = gene_sd[i] == 0
  )
}))

write.csv(
  long,
  "results/tables/exploratory_marker_expression.csv",
  row.names = FALSE
)

long$specimen_id <- factor(
  long$specimen_id, levels = colnames(norm)
)
long$gene_name <- factor(
  long$gene_name, levels = rev(genes)
)
long$group <- factor(long$group, levels = names(panels))

limit <- max(abs(long$row_z_score))
if (limit == 0) limit <- 1

p <- ggplot(
  long,
  aes(specimen_id, gene_name, fill = row_z_score)
) +
  geom_tile(colour = "white", linewidth = 0.3) +
  facet_grid(
    group ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    limits = c(-limit, limit),
    name = "Within-gene\nz-score"
  ) +
  labs(
    title = "Exploratory cell- and organ-associated gene expression",
    subtitle = "Paired specimens: Liver then Lung within each patient",
    x = NULL,
    y = NULL,
    caption = paste(
      "Row z-scores of log2(DESeq2-normalized count + 1).",
      "Colours compare specimens within each gene; they are not cell fractions.",
      "Marker expression does not establish cellular origin or tumor purity.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(size = 9),
    plot.caption = element_text(hjust = 0)
  )

ggsave(
  "results/figures/exploratory_marker_heatmap.png",
  p, width = 12, height = 10, dpi = 200
)
ggsave(
  "results/figures/exploratory_marker_heatmap.pdf",
  p, width = 12, height = 10
)

cat("\nMarker mapping checks passed; expression table and plots saved.\n")
