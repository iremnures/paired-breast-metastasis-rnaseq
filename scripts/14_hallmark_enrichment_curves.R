suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
})

analysis <- readRDS(
  "data/processed/hallmark_gsea_analysis.rds"
)

selected <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_COAGULATION",
  "HALLMARK_COMPLEMENT"
)

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

make_plot <- function(pathway, model) {
  result <- analysis$results[[model]]
  i <- match(pathway, result$pathway)
  stopifnot(!is.na(i))

  plotEnrichment(
    analysis$pathways[[pathway]],
    analysis$ranks[[model]],
    gseaParam = 1
  ) +
    labs(
      title = if (model == "primary") {
        "Primary: all six patients"
      } else {
        "Sensitivity: excluding R18"
      },
      subtitle = sprintf(
        "NES = %.2f | FDR = %.3g",
        result$NES[i], result$padj[i]
      ),
      x = "Gene rank: Lung direction (left) to Liver direction (right)",
      y = "Running enrichment score",
      caption = gsub(
        "_", " ",
        sub("^HALLMARK_", "", pathway)
      )
    ) +
    coord_cartesian(ylim = c(-1, 1)) +
    theme_bw(base_size = 11)
}

draw_pair <- function(pathway) {
  grid::grid.newpage()
  grid::pushViewport(
    grid::viewport(layout = grid::grid.layout(1, 2))
  )

  for (column in 1:2) {
    model <- c("primary", "no_R18")[column]
    print(
      make_plot(pathway, model),
      vp = grid::viewport(
        layout.pos.row = 1,
        layout.pos.col = column
      ),
      newpage = FALSE
    )
  }

  grid::popViewport()
}

pdf(
  "results/figures/hallmark_enrichment_curves.pdf",
  width = 13,
  height = 5
)
for (pathway in selected) {
  draw_pair(pathway)
}
dev.off()

for (pathway in selected) {
  png(
    paste0(
      "results/figures/enrichment_",
      sub("^HALLMARK_", "", pathway),
      ".png"
    ),
    width = 2600,
    height = 1000,
    res = 200
  )
  draw_pair(pathway)
  dev.off()
}

cat("\nEnrichment curves saved.\n")
