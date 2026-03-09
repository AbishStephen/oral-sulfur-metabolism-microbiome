source("scripts/01_load_data.R")
source("scripts/02_preprocess_data.R")

library(FactoMineR)
library(factoextra)
library(ggplot2)

cytokines <- merged_data[c(18:37)]

pca_results <- PCA(cytokines, scale.unit = TRUE)

# Order diagnosis groups
merged_data$diagnosis <- factor(merged_data$diagnosis, levels = c("LPH", "HPH", "G", "P"))

# Visualize PCA with specific requirements
pca_cytokine_plot <- fviz_pca_ind(pca_results, 
             geom.ind = "point", 
             col.ind = merged_data$diagnosis, 
             palette = c("gray", "#33FF57", "#3357FF", "#FF33A1"), # Bolder color palette
             addEllipses = TRUE, 
             legend.title = "Diagnosis",
             pointshape = 19,  # Use dots for all points
             pointsize = 3) +   # Fixed point size
  theme_minimal() +            # Use minimal theme for clean look
  theme(legend.position = "right", 
        legend.title = element_text(face = "bold"), 
        legend.text = element_text(size = 12), 
        panel.border = element_rect(color = "black", fill = NA, size = 1)) # Ensure color order
clinical_data$diagnosis <- factor(clinical_data$diagnosis,
                                  levels = c("LPH","HPH","G","P"))

clinical_data <- clinical_data %>%
  mutate(
    gcf_h2s_um_log = log1p(gcf_h2s_um),
    h2s_log = log1p(h2s),
    ch3sh_log = log1p(ch3sh)
  )

plot_gcf_h2s_um <- ggplot(clinical_data,
                          aes(x = diagnosis, y = gcf_h2s_um_log)) +
  geom_boxplot() +
  geom_jitter(width = 0.2)

ggsave("figures/figure1/sulcular_sulfide_diagnosis.svg",
       plot_gcf_h2s_um)
