library(DESeq2)

if (!requireNamespace(
  "ggplot2",
  quietly = TRUE
)) {
  install.packages("ggplot2")
}

if (!requireNamespace(
  "pheatmap",
  quietly = TRUE
)) {
  install.packages("pheatmap")
}

library(ggplot2)
library(pheatmap)

# ------------------------------------------------------------
# M07 - Sample-level QC
# Paired Liver vs Lung metastatic breast cancer
# ------------------------------------------------------------

# ------------------------------------------------------------
# 1. Read filtered biological-specimen counts
# ------------------------------------------------------------

count_df <- read.csv(
  "data/processed/primary_liver_lung_counts_filtered.csv.gz",
  check.names = FALSE
)

metadata <- read.csv(
  "data/processed/primary_liver_lung_metadata.csv",
  stringsAsFactors = FALSE
)

annotation_cols <- c(
  "gene_name",
  "Geneid",
  "Chr",
  "Start",
  "End",
  "Strand",
  "Length"
)

count_cols <- setdiff(
  colnames(count_df),
  annotation_cols
)

# ------------------------------------------------------------
# 2. Build count matrix
# ------------------------------------------------------------

count_matrix <- as.matrix(
  count_df[, count_cols]
)

rownames(count_matrix) <- count_df$Geneid

storage.mode(count_matrix) <- "integer"

# ------------------------------------------------------------
# 3. Reorder metadata exactly to count matrix
# ------------------------------------------------------------

metadata <- metadata[
  match(
    colnames(count_matrix),
    metadata$specimen_id
  ),
]

stopifnot(
  all(
    metadata$specimen_id
    == colnames(count_matrix)
  )
)

# ------------------------------------------------------------
# 4. Prepare design variables
# ------------------------------------------------------------

metadata$patient_code <- factor(
  metadata$patient_code
)

metadata$tissue <- factor(
  metadata$tissue,
  levels = c(
    "Liver",
    "Lung"
  )
)

rownames(metadata) <- metadata$specimen_id

# ------------------------------------------------------------
# 5. Create DESeq2 object
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ patient_code + tissue
)

# ------------------------------------------------------------
# 6. VST for QC / visualization only
# ------------------------------------------------------------

vsd <- vst(
  dds,
  blind = TRUE
)

vst_matrix <- assay(vsd)

# Save VST matrix
write.csv(
  vst_matrix,
  "data/processed/vst_matrix_qc.csv"
)

# ------------------------------------------------------------
# 7. Sample distance heatmap
# ------------------------------------------------------------

sample_dist <- dist(
  t(vst_matrix)
)

sample_dist_matrix <- as.matrix(
  sample_dist
)

annotation <- data.frame(
  Patient = metadata$patient_code,
  Tissue = metadata$tissue
)

rownames(annotation) <- metadata$specimen_id

png(
  "results/figures/01_sample_distance_heatmap.png",
  width = 1800,
  height = 1600,
  res = 200
)

pheatmap(
  sample_dist_matrix,
  annotation_col = annotation,
  annotation_row = annotation,
  main = "Sample-to-sample distances",
  border_color = NA
)

dev.off()

# ------------------------------------------------------------
# 8. PCA
# ------------------------------------------------------------

pca_data <- plotPCA(
  vsd,
  intgroup = c(
    "tissue",
    "patient_code"
  ),
  returnData = TRUE
)

percent_var <- round(
  100 * attr(
    pca_data,
    "percentVar"
  )
)

pca_data$specimen_id <- rownames(
  pca_data
)

write.csv(
  pca_data,
  "results/tables/pca_scores.csv",
  row.names = FALSE
)

p <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = tissue,
    label = specimen_id
  )
) +
  geom_point(
    size = 4
  ) +
  geom_text(
    vjust = -0.8,
    size = 3.4,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Liver" = "#7A1F3D",
      "Lung" = "#17365D"
    )
  ) +
  labs(
    title = "PCA of paired liver and lung metastases",
    x = paste0(
      "PC1: ",
      percent_var[1],
      "% variance"
    ),
    y = paste0(
      "PC2: ",
      percent_var[2],
      "% variance"
    ),
    color = "Tissue"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  "results/figures/02_pca_liver_lung.png",
  p,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 9. Print QC summary
# ------------------------------------------------------------

cat(
  "\nCount matrix dimensions:\n"
)

print(
  dim(count_matrix)
)

cat(
  "\nMetadata order matches count matrix:\n"
)

print(
  all(
    metadata$specimen_id
    == colnames(count_matrix)
  )
)

cat(
  "\nTissue counts:\n"
)

print(
  table(
    metadata$tissue
  )
)

cat(
  "\nPatients:\n"
)

print(
  table(
    metadata$patient_code
  )
)

cat(
  "\nPCA variance explained:\n"
)

print(
  percent_var
)

cat(
  "\nPCA scores:\n"
)

print(
  pca_data[
    ,
    c(
      "specimen_id",
      "patient_code",
      "tissue",
      "PC1",
      "PC2"
    )
  ]
)
