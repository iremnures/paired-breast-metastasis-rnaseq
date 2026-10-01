packages <- c(
  "DESeq2", "apeglm", "ggplot2", "ggrepel",
  "pheatmap", "fgsea", "msigdbr", "readxl",
  "BiocParallel"
)

available <- vapply(
  packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)

versions <- vapply(packages, function(package) {
  if (available[package]) {
    as.character(packageVersion(package))
  } else {
    NA_character_
  }
}, FUN.VALUE = character(1))

tab <- data.frame(
  package = packages,
  available = unname(available),
  version = unname(versions)
)

dir.create("docs", showWarnings = FALSE)

write.csv(
  tab,
  "docs/software_versions.csv",
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = "docs/environment_session_info.txt"
)

cat("\nR version:", R.version.string, "\n\n")
print(tab, row.names = FALSE)

if (any(!available)) {
  stop("Required packages are missing; inspect the table.")
}

cat("\nAll required packages are available.\n")
