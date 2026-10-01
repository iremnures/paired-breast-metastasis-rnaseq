suppressPackageStartupMessages(library(ggplot2))

tab <- read.csv(
  "results/tables/hallmark_primary_vs_no_R18.csv",
  check.names = FALSE
)

stopifnot(
  !anyDuplicated(tab$pathway),
  all(is.finite(tab$primary_NES)),
  all(is.finite(tab$no_R18_NES))
)

# Sort by primary NES: Liver-direction at the bottom,
# Lung-direction at the top.
tab <- tab[order(tab$primary_NES, tab$pathway), ]

labels <- tools::toTitleCase(
  tolower(gsub("_", " ", sub("^HALLMARK_", "", tab$pathway)))
)

label_levels <- labels

tab$pathway_label <- factor(
  labels,
  levels = label_levels
)

long <- rbind(
  data.frame(
    pathway_label = tab$pathway_label,
    model = "Primary",
    NES = tab$primary_NES,
    padj = tab$primary_padj
  ),
  data.frame(
    pathway_label = tab$pathway_label,
    model = "Excluding R18",
    NES = tab$no_R18_NES,
    padj = tab$no_R18_padj
  )
)

long$pathway_label <- factor(
  long$pathway_label,
  levels = label_levels
)

long$model <- factor(
  long$model,
  levels = c("Primary", "Excluding R18")
)

long$FDR_significant <- (
  is.finite(long$padj) & long$padj < 0.05
)

p <- ggplot() +
  geom_segment(
    data = tab,
    aes(
      x = primary_NES,
      xend = no_R18_NES,
      y = pathway_label,
      yend = pathway_label
    ),
    colour = "#B8B8B8",
    linewidth = 0.6
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "#666666"
  ) +
  geom_point(
    data = long,
    aes(
      x = NES,
      y = pathway_label,
      shape = model,
      colour = model,
      fill = model
    ),
    size = 2.8,
    stroke = 0.8
  ) +
  # White centres identify results without FDR < 0.05 support.
  geom_point(
    data = long[!long$FDR_significant, ],
    aes(
      x = NES,
      y = pathway_label,
      shape = model,
      colour = model
    ),
    fill = "white",
    size = 2.8,
    stroke = 0.8,
    show.legend = FALSE
  ) +
  scale_shape_manual(
    values = c("Primary" = 21, "Excluding R18" = 24)
  ) +
  scale_colour_manual(
    values = c("Primary" = "#007C83", "Excluding R18" = "#D87521")
  ) +
  scale_fill_manual(
    values = c("Primary" = "#007C83", "Excluding R18" = "#D87521")
  ) +
  labs(
    title = "Hallmark enrichment: primary and R18 sensitivity models",
    subtitle = "Identical gene universe and gene-set membership in both models",
    x = "Normalized enrichment score (NES)\nLiver direction < 0 | Lung direction > 0",
    y = NULL,
    shape = "Model",
    colour = "Model",
    fill = "Model",
    caption = paste(
      "Filled symbols: FDR < 0.05. Open symbols: FDR >= 0.05 or unavailable.",
      "Connecting lines show NES differences, not confidence intervals.",
      "BH correction was applied separately within each model.",
      sep = "\n"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_line(colour = "#EEEEEE"),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.y = element_text(size = 9),
    plot.caption = element_text(hjust = 0)
  )

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "results/figures/hallmark_primary_vs_no_R18_summary.png",
  plot = p,
  width = 12,
  height = 15,
  dpi = 200
)

ggsave(
  "results/figures/hallmark_primary_vs_no_R18_summary.pdf",
  plot = p,
  width = 12,
  height = 15
)

cat("\nHallmark summary plot saved.\n")
