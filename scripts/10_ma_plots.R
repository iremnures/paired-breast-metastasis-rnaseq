suppressPackageStartupMessages(library(DESeq2))

tab <- read.csv(
  "results/tables/apeglm_results_with_diagnostic_flags.csv",
  check.names = FALSE
)

dds <- readRDS("data/processed/deseq2_primary_dds.rds")

raw <- results(
  dds,
  contrast = c("tissue", "Lung", "Liver"),
  alpha = 0.05
)

idx <- match(tab$Geneid, rownames(raw))
stopifnot(!anyNA(idx))

x <- log10(raw$baseMean[idx] + 1)
sig <- !is.na(tab$primary_padj) & tab$primary_padj < 0.05
flag <- tab$optimizer_sensitive

colors <- rep("#BDBDBD", nrow(tab))
colors[sig] <- "#007C83"
colors[flag] <- "#D87521"

y_values <- c(tab$log2FC_MLE, tab$default_apeglm)
limit <- max(abs(y_values[is.finite(y_values)])) * 1.05

draw <- function() {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

  for (column in c("log2FC_MLE", "default_apeglm")) {
    y <- tab[[column]]
    valid <- is.finite(x) & is.finite(y)

    # Draw flagged genes last so they remain visible.
    order_idx <- c(
      which(valid & !flag),
      which(valid & flag)
    )

    plot(
      x[order_idx], y[order_idx],
      pch = 16, cex = 0.45,
      col = colors[order_idx],
      ylim = c(-limit, limit),
      xlab = "log10(mean normalized count + 1)",
      ylab = "log2 fold change: Lung vs Liver",
      main = if (column == "log2FC_MLE") {
        "Before shrinkage"
      } else {
        "Default apeglm: diagnostic view"
      },
      bty = "l"
    )

    abline(h = 0, col = "#555555", lty = 2)

    legend(
      "topleft",
      legend = c(
        "Primary FDR < 0.05",
        "Other genes",
        "Optimizer-sensitive"
      ),
      col = c("#007C83", "#BDBDBD", "#D87521"),
      pch = 16, bty = "n", cex = 0.8
    )
  }
}

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

pdf("results/figures/ma_before_after_shrinkage.pdf",
    width = 12, height = 6)
draw()
dev.off()

png("results/figures/ma_before_after_shrinkage.png",
    width = 2400, height = 1200, res = 200)
draw()
dev.off()

cat("\nMA plots saved.\n")
