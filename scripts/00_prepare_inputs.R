suppressPackageStartupMessages(library(readxl))

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) {
  args[1]
} else {
  "data/processed/rebuild_check"
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read.csv(
  "data/raw/GSE316391_counts_PE.csv.gz",
  check.names = FALSE
)

rna <- as.data.frame(read_excel(
  "data/metadata/GSE316391_clinical_metadata.xlsx",
  sheet = "RNA-seq",
  col_types = "text"
))

annotation_cols <- c(
  "gene_name", "Geneid", "Chr",
  "Start", "End", "Strand", "Length"
)

required <- c(
  "RNA-seq name", "Patient #", "Tissue Site",
  "External ID", "Container ID", "Internal Case ID",
  "Histology"
)

stopifnot(
  all(annotation_cols %in% names(raw)),
  all(required %in% names(rna)),
  !anyNA(raw$Geneid),
  !anyDuplicated(raw$Geneid)
)

library_cols <- setdiff(names(raw), annotation_cols)

stopifnot(
  !anyDuplicated(library_cols),
  all(grepl("^R[0-9]+(BR|LIV|LUN)_[0-9]+$", library_cols))
)

# Normalize only the final library-replicate separator.
rna$library_name <- sub(
  "[-_]([0-9]+)$", "_\\1",
  trimws(rna[["RNA-seq name"]])
)

rna <- rna[
  !is.na(rna$library_name) &
  rna$library_name %in% library_cols,
  , drop = FALSE
]

stopifnot(
  !anyDuplicated(rna$library_name),
  setequal(rna$library_name, library_cols)
)

rna$specimen_id <- sub("_[0-9]+$", "", rna$library_name)

count_matrix <- as.matrix(raw[, library_cols, drop = FALSE])

stopifnot(
  is.numeric(count_matrix),
  all(is.finite(count_matrix)),
  all(count_matrix >= 0),
  all(count_matrix == floor(count_matrix))
)

# Preserve specimen order from the raw-count columns.
count_specimens <- unique(sub("_[0-9]+$", "", library_cols))

summed <- vapply(
  count_specimens,
  function(specimen) {
    libraries <- library_cols[
      sub("_[0-9]+$", "", library_cols) == specimen
    ]
    rowSums(count_matrix[, libraries, drop = FALSE])
  },
  FUN.VALUE = numeric(nrow(raw))
)

colnames(summed) <- count_specimens
stopifnot(all(summed <= .Machine$integer.max))

all_counts <- cbind(
  raw[, annotation_cols, drop = FALSE],
  as.data.frame(summed, check.names = FALSE)
)

# Metadata order follows the source RNA-seq worksheet.
metadata_specimens <- unique(rna$specimen_id)
metadata_rows <- list()

for (specimen in metadata_specimens) {
  rows <- rna[rna$specimen_id == specimen, , drop = FALSE]

  # Library records for one specimen must agree on its metadata.
  for (column in setdiff(required, "RNA-seq name")) {
    values <- trimws(rows[[column]])
    stopifnot(length(unique(values)) == 1L)
  }

  patient <- sub("(BR|LIV|LUN)$", "", specimen)
  tissue_code <- sub("^R[0-9]+", "", specimen)

  tissue <- unname(c(
    BR = "Breast", LIV = "Liver", LUN = "Lung"
  )[tissue_code])

  clinical_site <- trimws(rows[["Tissue Site"]][1])
  expected_site <- unname(c(
    BR = "Breast", LIV = "Liver/Bile Duct", LUN = "Lung"
  )[tissue_code])

  stopifnot(identical(clinical_site, expected_site))

  metadata_rows[[specimen]] <- data.frame(
    specimen_id = specimen,
    patient_code = patient,
    patient_number = as.integer(rows[["Patient #"]][1]),
    tissue_code = tissue_code,
    tissue = tissue,
    external_id = trimws(rows[["External ID"]][1]),
    container_id = trimws(rows[["Container ID"]][1]),
    internal_case_id = trimws(rows[["Internal Case ID"]][1]),
    histology = trimws(rows[["Histology"]][1]),
    n_rna_libraries = nrow(rows),
    library_names = paste(sort(rows$library_name), collapse = ",")
  )
}

metadata <- do.call(rbind, metadata_rows)
rownames(metadata) <- NULL

metastasis_meta <- metadata[
  metadata$tissue %in% c("Liver", "Lung"),
  , drop = FALSE
]

pairs <- table(
  metastasis_meta$patient_code,
  factor(metastasis_meta$tissue, levels = c("Liver", "Lung"))
)

paired_patients <- rownames(pairs)[
  pairs[, "Liver"] == 1 & pairs[, "Lung"] == 1
]

primary_meta <- metastasis_meta[
  metastasis_meta$patient_code %in% paired_patients,
  , drop = FALSE
]

primary_meta <- primary_meta[
  order(
    primary_meta$patient_code,
    match(primary_meta$tissue, c("Liver", "Lung"))
  ),
  , drop = FALSE
]
rownames(primary_meta) <- NULL

stopifnot(
  length(paired_patients) == 6,
  nrow(primary_meta) == 12
)

primary_counts <- all_counts[
  , c(annotation_cols, primary_meta$specimen_id),
  drop = FALSE
]

keep <- rowSums(
  as.matrix(primary_counts[, primary_meta$specimen_id]) >= 10
) >= 6

filtered_counts <- primary_counts[keep, , drop = FALSE]

outputs <- list(
  "counts_by_biological_specimen.csv.gz" = all_counts,
  "specimen_metadata.csv" = metadata,
  "primary_liver_lung_counts.csv.gz" = primary_counts,
  "primary_liver_lung_metadata.csv" = primary_meta,
  "primary_liver_lung_counts_filtered.csv.gz" = filtered_counts
)

write_table <- function(tab, path) {
  connection <- if (grepl("\\.gz$", path)) {
    gzfile(path, "wt")
  } else {
    file(path, "wt")
  }
  on.exit(close(connection))
  write.csv(tab, connection, row.names = FALSE)
}

checks <- list()

for (filename in names(outputs)) {
  rebuilt <- outputs[[filename]]
  write_table(rebuilt, file.path(output_dir, filename))

  existing_path <- file.path("data/processed", filename)
  matches <- NA

  if (file.exists(existing_path)) {
    existing <- read.csv(
      existing_path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    # Compare CSV content after writing and reading.
    rebuilt_read <- read.csv(
      file.path(output_dir, filename),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    matches <- isTRUE(all.equal(
      existing,
      rebuilt_read,
      check.attributes = FALSE,
      tolerance = 0
    ))

    if (!matches) {
      cat("\n--- Difference:", filename, "---\n")
      print(all.equal(
        existing,
        rebuilt_read,
        check.attributes = FALSE,
        tolerance = 0
      ))
    }
  }

  checks[[filename]] <- data.frame(
    file = filename,
    rows = nrow(rebuilt),
    columns = ncol(rebuilt),
    matches_existing = matches
  )
}

check_table <- do.call(rbind, checks)
rownames(check_table) <- NULL

write.csv(
  check_table,
  file.path(output_dir, "rebuild_comparison.csv"),
  row.names = FALSE
)

cat("\n--- Rebuilt inputs ---\n")
print(check_table, row.names = FALSE)

cat("\nPaired patients:\n")
print(paired_patients)

cat("\nGenes before filtering:", nrow(primary_counts), "\n")
cat("Genes retained:", nrow(filtered_counts), "\n")
cat("Output directory:", output_dir, "\n")

if (any(check_table$matches_existing %in% FALSE)) {
  stop("Rebuilt inputs differ from existing inputs; inspect before replacing.")
}
