suppressPackageStartupMessages(library(ggplot2))

analysis <- readRDS(
  "data/processed/leave_one_out/hallmark_analysis.rds"
)

tab <- analysis$summary
tab <- tab[order(tab$primary_NES, tab$pathway), ]

models <- c(
  "primary", "no_R18", "no_R23",
  "no_R33", "no_R36", "no_R49", "no_R8"
)

labels <- gsub("_", " ", sub("^HALLMARK_", "", tab$pathway))
long_list <- list()

for (model in models) {
  result <- analysis$results[[model]]
  idx <- match(tab$pathway, result$pathway)
  stopifnot(!anyNA(idx))

  long_list[[model]] <- data.frame(
    pathway = tab$pathway,
    pathway_label = labels,
    model = model,
    NES = result$NES[idx],
    padj = result$padj[idx]
  )
}

long <- do.call(rbind, long_list)
rownames(long) <- NULL

long$pathway_label <- factor(
  long$pathway_label,
  levels = labels
)
long$model <- factor(long$model, levels = models)
long$significant <- is.finite(long$padj) & long$padj < 0.05

limit <- max(abs(long$NES), na.rm = TRUE)

p <- ggplot(
  long,
  aes(x = model, y = pathway_label, fill = NES)
) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(
    data = long[long$significant, ],
    shape = 21,
    fill = "black",
    colour = "white",
    size = 1.8,
    stroke = 0.25
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-limit, limit),
    name = "NES"
  ) +
  scale_x_discrete(
    labels = c(
      "primary" = "Primary",
      "no_R18" = "No R18",
      "no_R23" = "No R23",
      "no_R33" = "No R33",
      "no_R36" = "No R36",
      "no_R49" = "No R49",
      "no_R8" = "No R8"
    )
  ) +
  labs(
    title = "Hallmark enrichment: leave-one-patient-out sensitivity",
    subtitle = "Blue: Liver direction | Red: Lung direction",
    x = NULL,
    y = NULL,
    caption = paste(
      "Dot: FDR < 0.05 within that model. No dot: threshold not met.",
      "All models use the same 25,667 ranked genes and 50 gene sets.",
      "These overlapping subsets are sensitivity analyses, not independent validations.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.caption = element_text(hjust = 0)
  )

ggsave(
  "results/figures/hallmark_leave_one_out_heatmap.png",
  p, width = 12, height = 15, dpi = 200
)
ggsave(
  "results/figures/hallmark_leave_one_out_heatmap.pdf",
  p, width = 12, height = 15
)

write.csv(
  long,
  "results/tables/leave_one_out/hallmark_heatmap_values.csv",
  row.names = FALSE
)

# Report every primary-significant pathway/model combination
# where direction changes or the FDR threshold is not retained.
primary_idx <- match(long$pathway, tab$pathway)
primary_sig <- tab$primary_FDR_below_0.05[primary_idx]

long$same_direction <- (
  sign(long$NES) == sign(tab$primary_NES[primary_idx])
)

issues <- long[
  long$model != "primary" &
  primary_sig &
  (!long$significant | !long$same_direction),
  c("pathway", "model", "NES", "padj", "same_direction")
]

write.csv(
  issues,
  "results/tables/leave_one_out/hallmark_sensitivity_changes.csv",
  row.names = FALSE
)

cat("\n--- Primary-significant pathways: sensitivity changes ---\n")
print(issues, row.names = FALSE)
cat("\nHeatmap saved.\n")
