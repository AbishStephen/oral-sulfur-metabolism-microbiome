source("scripts/01_load_data.R")
source("scripts/02_preprocess_data.R")

# ================================
# Figure 1
# ================================

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggplot2)

# ----------------
# Load data
# ----------------

clinical_data <- read_csv("data/clinical_data.csv")
cytokine_data <- read_csv("data/gcf_cytokines.csv")

# ----------------
# Merge data for PCA
# ----------------

clinical_data <- clinical_data %>%
  rownames_to_column(var = "patient_id")

cytokine_data <- cytokine_data %>%
  rownames_to_column(var = "patient_id")

merged_data <- clinical_data %>%
  inner_join(cytokine_data, by = "patient_id")

# ----------------
# PCA on cytokines
# ----------------

cytokines <- merged_data[c(18:37)]

pca_results <- PCA(cytokines, scale.unit = TRUE)

# ----------------
# PCA loading plot
# ----------------

cytokine_loadings <- pca_results$var$cos2
cytokine_loadings <- as.data.frame(cytokine_loadings)

cytokine_loadings$max_cos2 <- apply(
  cytokine_loadings[, c("Dim.1", "Dim.2")],
  1,
  max
)

p <- ggplot(cytokine_loadings,
            aes(x = Dim.1, y = Dim.2, color = max_cos2)) +
  
  geom_segment(
    aes(x = 0, xend = Dim.1, y = 0, yend = Dim.2),
    arrow = arrow(type = "open", length = unit(0.05, "inches")),
    size = 0.6
  ) +
  
  geom_text(aes(label = rownames(cytokine_loadings)),
            vjust = 1.5) +
  
  scale_color_gradient2(
    low = "darkgrey",
    mid = "orange",
    high = "black",
    midpoint = 0.02
  ) +
  
  theme_minimal() +
  
  theme(
    aspect.ratio = 1,
    axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "PC1",
    y = "PC2",
    color = "Cos2 Value"
  ) +
  
  coord_fixed(ratio = 1)

# ----------------
# PCA score plot
# ----------------

merged_data$diagnosis <- factor(
  merged_data$diagnosis,
  levels = c("LPH", "HPH", "G", "P")
)

pca_cytokine_plot <- fviz_pca_ind(
  pca_results,
  geom.ind = "point",
  col.ind = merged_data$diagnosis,
  palette = c("gray", "#33FF57", "#3357FF", "#FF33A1"),
  addEllipses = TRUE,
  legend.title = "Diagnosis",
  pointshape = 19,
  pointsize = 3
) +
  
  theme_minimal() +
  
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

ggsave(
  "figures/figure1/cytokine_pca_plot.svg",
  plot = pca_cytokine_plot,
  width = 10,
  height = 10,
  device = "svg"
)

# ===============================
# VSC boxplots
# ===============================

# Ensure the diagnosis column is a factor and ordered correctly
clinical_data$diagnosis <- factor(
  clinical_data$diagnosis,
  levels = c("LPH", "HPH", "G", "P")
)

# Filter out rows with NA in the diagnosis column
clinical_data <- clinical_data %>%
  filter(!is.na(diagnosis))

# Apply log transformation
clinical_data <- clinical_data %>%
  mutate(
    gcf_h2s_um_log = log1p(gcf_h2s_um),
    h2s_log = log1p(h2s),
    ch3sh_log = log1p(ch3sh)
  )

# ----------------
# Sulcular H2S
# ----------------

plot_gcf_h2s_um <- ggplot(
  clinical_data,
  aes(x = diagnosis, y = gcf_h2s_um_log, color = diagnosis)
) +
  
  geom_boxplot(outlier.shape = NA) +
  
  geom_jitter(width = 0.2, size = 2) +
  
  theme_minimal() +
  
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  
  labs(
    x = "Diagnosis",
    y = expression(Log ~ "[" ~ H[2] ~ S ~ "]" ~ (mu ~ M))
  ) +
  
  scale_color_manual(
    values = c("gray", "#33FF57", "#3357FF", "#FF33A1")
  )

# ----------------
# Breath H2S
# ----------------

plot_h2s <- ggplot(
  clinical_data,
  aes(x = diagnosis, y = h2s_log, color = diagnosis)
) +
  
  geom_boxplot(outlier.shape = NA) +
  
  geom_jitter(width = 0.2, size = 2) +
  
  theme_minimal() +
  
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  
  labs(
    x = "Diagnosis",
    y = expression(Log ~ "[" ~ H[2] ~ S ~ "]" ~ (ppb))
  ) +
  
  scale_color_manual(
    values = c("gray", "#33FF57", "#3357FF", "#FF33A1")
  )

# ----------------
# Breath CH3SH
# ----------------

plot_ch3sh <- ggplot(
  clinical_data,
  aes(x = diagnosis, y = ch3sh_log, color = diagnosis)
) +
  
  geom_boxplot(outlier.shape = NA) +
  
  geom_jitter(width = 0.2, size = 2) +
  
  theme_minimal() +
  
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  
  labs(
    x = "Diagnosis",
    y = expression(Log ~ "[" ~ CH[3] ~ SH ~ "]" ~ (ppb))
  ) +
  
  scale_color_manual(
    values = c("gray", "#33FF57", "#3357FF", "#FF33A1")
  )

# ----------------
# Print plots
# ----------------

print(plot_gcf_h2s_um)
print(plot_h2s)
print(plot_ch3sh)

# ----------------
# Save plots
# ----------------

ggsave(
  "figures/figure1/sulcular_sulfide_diagnosis.svg",
  plot = plot_gcf_h2s_um,
  width = 5,
  height = 5,
  device = "svg"
)

ggsave(
  "figures/figure1/breath_h2s_diagnosis.svg",
  plot = plot_h2s,
  width = 5,
  height = 5,
  device = "svg"
)

ggsave(
  "figures/figure1/breath_ch3sh_diagnosis.svg",
  plot = plot_ch3sh,
  width = 5,
  height = 5,
  device = "svg"
)
