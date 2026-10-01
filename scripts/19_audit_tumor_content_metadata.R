suppressPackageStartupMessages({
  library(DESeq2)
  library(readxl)
})

dds <- readRDS(
  "data/processed/deseq2_primary_dds.rds"
)
meta <- as.data.frame(colData(dds))

clinical <- as.data.frame(read_excel(
  "data/metadata/GSE316391_clinical_metadata.xlsx",
  sheet = "Clinical Data",
  col_types = "text",
  .name_repair = "unique"
))

required <- c(
  "Study ID",
  "Container ID",
  "Internal Case ID",
  "Tissue Site",
  "Primary Diagnosis",
  "Percent Tumor In Biopsy",
  "Percent Necrosis In Biopsy"
)

# Verify and select the first specimen-identity block.
stopifnot(
  names(clinical)[2] == "Tissue Site...2",
  names(clinical)[3] == "Study ID...3"
)

names(clinical)[2:3] <- c("Tissue Site", "Study ID")

stopifnot(all(required %in% names(clinical)))

clean <- function(x) trimws(as.character(x))

for (column in required) {
  clinical[[column]] <- clean(clinical[[column]])
}

# Ignore blank separator rows.
clinical <- clinical[
  !is.na(clinical[["Container ID"]]) &
  nzchar(clinical[["Container ID"]]),
  , drop = FALSE
]

key <- function(study, container, case) {
  paste(clean(study), clean(container), clean(case), sep = "|")
}

clinical_key <- key(
  clinical[["Study ID"]],
  clinical[["Container ID"]],
  clinical[["Internal Case ID"]]
)

metadata_key <- key(
  meta$external_id,
  meta$container_id,
  meta$internal_case_id
)

stopifnot(
  !anyDuplicated(clinical_key),
  !anyDuplicated(metadata_key)
)

idx <- match(metadata_key, clinical_key)

if (anyNA(idx)) {
  print(meta[is.na(idx), c(
    "specimen_id", "external_id",
    "container_id", "internal_case_id"
  )])
  stop("Some specimens could not be matched; inspect before proceeding.")
}

audit <- data.frame(
  specimen_id = meta$specimen_id,
  patient_code = meta$patient_code,
  tissue = meta$tissue,
  clinical_tissue = clinical[["Tissue Site"]][idx],
  clinical_histology = clinical[["Primary Diagnosis"]][idx],
  tumor_content_original = clinical[["Percent Tumor In Biopsy"]][idx],
  necrosis_original = clinical[["Percent Necrosis In Biopsy"]][idx]
)

missing_record <- function(x) {
  is.na(x) |
    tolower(trimws(x)) %in% c(
      "", "unknown", "n/a", "na", "not available", "not reported"
    )
}

audit$tumor_content_missing <- missing_record(
  audit$tumor_content_original
)
audit$necrosis_missing <- missing_record(
  audit$necrosis_original
)

# Preserve ranges and other recorded descriptions without conversion.
write.csv(
  audit,
  "results/tables/tumor_content_metadata_audit.csv",
  row.names = FALSE
)

cat("\n--- Matched clinical records: all 12 specimens ---\n")
print(audit[, c(
  "specimen_id", "tissue", "clinical_tissue",
  "tumor_content_original", "necrosis_original"
)], row.names = FALSE)

cat("\n--- Tumor-content records ---\n")
print(table(
  audit$tumor_content_original,
  useNA = "ifany"
))

cat("\n--- Missing tumor-content records by site ---\n")
print(with(
  audit,
  table(tissue, tumor_content_missing)
))

cat("\nSpecimens with a recorded tumor-content description:",
    sum(!audit$tumor_content_missing), "of", nrow(audit), "\n")

cat("\nPairs with tumor-content descriptions at both sites:\n")
pair_complete <- tapply(
  !audit$tumor_content_missing,
  audit$patient_code,
  all
)
print(pair_complete)

cat("\n--- Histology of matched specimens ---\n")
print(audit[, c(
  "specimen_id", "clinical_histology"
)], row.names = FALSE)
