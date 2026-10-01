library(DESeq2)

# ============================================================
# M07 - Paired DESeq2 model
# Liver vs Lung metastases
# ============================================================

# ------------------------------------------------------------
# 1. Read filtered count matrix and metadata
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

sample_cols <- setdiff(
  colnames(count_df),
  annotation_cols
)

# ------------------------------------------------------------
# 2. Build raw count matrix
# ------------------------------------------------------------

count_matrix <- as.matrix(
  count_df[, sample_cols]
)

rownames(count_matrix) <- count_df$Geneid

storage.mode(count_matrix) <- "integer"

# ------------------------------------------------------------
# 3. Match metadata to counts
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

rownames(metadata) <- metadata$specimen_id

# ------------------------------------------------------------
# 4. Define factors explicitly
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

# ------------------------------------------------------------
# 5. Check paired structure
# ------------------------------------------------------------

cat("\nPatient x Tissue table:\n")

print(
  table(
    metadata$patient_code,
    metadata$tissue
  )
)

# ------------------------------------------------------------
# 6. Inspect design matrix
# ------------------------------------------------------------

design_matrix <- model.matrix(
  ~ patient_code + tissue,
  data = metadata
)

cat("\nDesign matrix:\n")
print(design_matrix)

cat("\nDesign matrix rank:\n")
print(qr(design_matrix)$rank)

cat("\nNumber of design columns:\n")
print(ncol(design_matrix))

cat("\nFull rank:\n")
print(
  qr(design_matrix)$rank
  == ncol(design_matrix)
)

# ------------------------------------------------------------
# 7. Create DESeq2 dataset
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ patient_code + tissue
)

# ------------------------------------------------------------
# 8. Run DESeq2
# ------------------------------------------------------------

dds <- DESeq(
  dds
)

# ------------------------------------------------------------
# 9. Inspect coefficient names
# ------------------------------------------------------------

cat("\nDESeq2 coefficient names:\n")
print(
  resultsNames(dds)
)

# ------------------------------------------------------------
# 10. Size factors
# ------------------------------------------------------------

size_factor_table <- data.frame(
  specimen_id = colnames(dds),
  patient = colData(dds)$patient_code,
  tissue = colData(dds)$tissue,
  raw_library_size = colSums(
    counts(dds)
  ),
  size_factor = sizeFactors(dds)
)

cat("\nSize factors:\n")

print(
  size_factor_table
)

write.csv(
  size_factor_table,
  "results/tables/deseq2_size_factors.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Save normalized counts
# ------------------------------------------------------------

normalized_counts <- counts(
  dds,
  normalized = TRUE
)

write.csv(
  normalized_counts,
  "data/processed/deseq2_normalized_counts.csv"
)

# ------------------------------------------------------------
# 12. Save fitted DESeq2 object
# ------------------------------------------------------------

saveRDS(
  dds,
  "data/processed/deseq2_primary_dds.rds"
)

# ------------------------------------------------------------
# 13. Session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = "docs/sessionInfo_M07.txt"
)

cat("\nDESeq2 model fitting completed successfully.\n")
