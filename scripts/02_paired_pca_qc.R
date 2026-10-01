library(ggplot2)
library(ggrepel)

pca <- read.csv(
  "results/tables/pca_scores.csv",
  stringsAsFactors = FALSE
)

# Ensure paired order is explicit
pca$tissue <- factor(
  pca$tissue,
  levels = c("Liver", "Lung")
)

# ------------------------------------------------------------
# Paired PCA
# ------------------------------------------------------------

p <- ggplot(
  pca,
  aes(
    x = PC1,
    y = PC2
  )
) +

  # Connect matched specimens from the same patient
  geom_line(
    aes(
      group = patient_code
    ),
    linewidth = 0.7,
    alpha = 0.55
  ) +

  geom_point(
    aes(
      color = tissue
    ),
    size = 4
  ) +

  geom_text_repel(
    aes(
      label = specimen_id,
      color = tissue
    ),
    size = 3.5,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +

  scale_color_manual(
    values = c(
      "Liver" = "#7A1F3D",
      "Lung" = "#17365D"
    )
  ) +

  labs(
    title = "Paired PCA of liver and lung metastases",
    subtitle = "Lines connect metastatic sites from the same patient",
    x = "PC1: 37% variance",
    y = "PC2: 18% variance",
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
  "results/figures/03_paired_pca_liver_lung.png",
  p,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

# ------------------------------------------------------------
# Paired PCA differences
# ------------------------------------------------------------

wide_pc1 <- reshape(
  pca[
    ,
    c(
      "patient_code",
      "tissue",
      "PC1"
    )
  ],
  idvar = "patient_code",
  timevar = "tissue",
  direction = "wide"
)

wide_pc2 <- reshape(
  pca[
    ,
    c(
      "patient_code",
      "tissue",
      "PC2"
    )
  ],
  idvar = "patient_code",
  timevar = "tissue",
  direction = "wide"
)

paired <- merge(
  wide_pc1,
  wide_pc2,
  by = "patient_code"
)

paired$PC1_Lung_minus_Liver <- (
  paired$PC1.Lung -
  paired$PC1.Liver
)

paired$PC2_Lung_minus_Liver <- (
  paired$PC2.Lung -
  paired$PC2.Liver
)

print(
  paired[
    ,
    c(
      "patient_code",
      "PC1_Lung_minus_Liver",
      "PC2_Lung_minus_Liver"
    )
  ]
)

write.csv(
  paired,
  "results/tables/paired_pca_differences.csv",
  row.names = FALSE
)
