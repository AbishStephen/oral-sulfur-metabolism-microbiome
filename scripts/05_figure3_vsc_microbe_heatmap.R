# ========================================
# Figure 3
# Microbial taxa correlated with VSCs
# ========================================

library(tidyverse)
library(reshape2)
library(ggplot2)
library(psych)
library(patchwork)
library(compositions)
library(pheatmap)
library(RColorBrewer)
library(ggh4x)
library(ggnewscale)

# ------------------------------
# Load data
# ------------------------------

sample_info <- read.csv("data/sampleinfo.csv")
clinical_data <- read.csv("data/clinical_data.csv")
gcf_cytokines <- read.csv("data/gcf_cytokines.csv")
microbial_data <- read.csv("data/microbial_data.csv")

# Remove duplicated diagnosis column
sample_info <- sample_info[-c(4)]

# ------------------------------
# Merge datasets
# ------------------------------

microbial_with_info <- microbial_data %>%
  left_join(sample_info, by = "sample")

merged_data <- microbial_with_info %>%
  left_join(clinical_data, by = "patno") %>%
  left_join(gcf_cytokines, by = "patno")

# ------------------------------
# Transform variables
# ------------------------------

merged_data <- merged_data %>%
  mutate(
    h2s_log = log1p(h2s),
    ch3sh_log = log1p(ch3sh),
    gcf_h2s_um_log = log1p(gcf_h2s_um),
    gcf_ch3sh_um_log = log1p(gcf_ch3sh_um)
  )

cytokine_cols <- names(gcf_cytokines)[!(names(gcf_cytokines) %in% c("patno"))]

merged_data[cytokine_cols] <- scale(merged_data[cytokine_cols])

microbe_cols <- setdiff(names(microbial_data), "sample")

merged_data[microbe_cols] <- merged_data[microbe_cols] + 1
merged_data[microbe_cols] <- clr(merged_data[microbe_cols])

# ------------------------------
# Correlation analysis
# ------------------------------

niches <- c("interdental", "subgingival", "tongue")

significant_correlations <- list()

for (niche in niches) {
  
  niche_data <- merged_data %>%
    filter(niche == !!niche) %>%
    filter(!is.na(diagnosis))
  
  diagnosis_groups <- unique(niche_data$diagnosis)
  
  for (diagnosis in diagnosis_groups) {
    
    diag_data <- niche_data %>%
      filter(diagnosis == !!diagnosis)
    
    selected_data <- diag_data %>%
      select(
        h2s_log,
        ch3sh_log,
        gcf_h2s_um_log,
        gcf_ch3sh_um_log,
        all_of(cytokine_cols),
        all_of(microbe_cols)
      )
    
    if (nrow(selected_data) > 5) {
      
      cor_matrix <- corr.test(selected_data, method = "spearman")
      
      cor_coeffs <- cor_matrix$r
      p_values <- cor_matrix$p
      
      p_adj <- p.adjust(p_values, method = "fdr")
      
      cor_coeffs_filtered <- cor_coeffs
      
      cor_coeffs_filtered[p_adj > 0.05 | abs(cor_coeffs_filtered) < 0.4] <- NA
      
      sig_pairs <- which(!is.na(cor_coeffs_filtered), arr.ind = TRUE)
      
      sig_correlations <- data.frame(
        from = rownames(cor_coeffs_filtered)[sig_pairs[, 1]],
        to = colnames(cor_coeffs_filtered)[sig_pairs[, 2]],
        weight = cor_coeffs_filtered[sig_pairs]
      )
      
      sig_correlations <- sig_correlations %>%
        filter(from != to) %>%
        filter(!(
          (from %in% microbe_cols & to %in% microbe_cols) |
          (from %in% cytokine_cols & to %in% cytokine_cols) |
          (from %in% c("h2s_log","ch3sh_log","gcf_h2s_um_log","gcf_ch3sh_um_log") &
           to %in% c("h2s_log","ch3sh_log","gcf_h2s_um_log","gcf_ch3sh_um_log"))
        )) %>%
        mutate(
          niche = niche,
          diagnosis = diagnosis
        )
      
      significant_correlations[[paste(niche, diagnosis, sep = "_")]] <- sig_correlations
      
    }
  }
}

# ------------------------------
# Extract microbes correlated with VSCs
# ------------------------------

vsc_targets <- c("gcf_h2s_um_log","gcf_ch3sh_um_log")

microbes_vsc_by_niche <- list()

for (niche in niches) {
  
  niche_keys <- names(significant_correlations)[grepl(paste0("^", niche, "_"), names(significant_correlations))]
  
  microbes_in_niche <- c()
  
  for (key in niche_keys) {
    
    sig_corrs <- significant_correlations[[key]]
    
    if (!is.null(sig_corrs) && nrow(sig_corrs) > 0) {
      
      vsc_microbe_corrs <- sig_corrs %>%
        filter(
          (from %in% microbe_cols & to %in% vsc_targets) |
          (to %in% microbe_cols & from %in% vsc_targets)
        )
      
      microbes_found <- vsc_microbe_corrs %>%
        mutate(microbe = ifelse(from %in% microbe_cols, from, to)) %>%
        pull(microbe)
      
      microbes_in_niche <- c(microbes_in_niche, microbes_found)
      
    }
  }
  
  microbes_vsc_by_niche[[niche]] <- unique(microbes_in_niche)
}

# ------------------------------
# Prepare CLR abundance heatmap
# ------------------------------

vsc_microbe_combined <- bind_rows(significant_correlations, .id = "source") %>%
  filter(
    (from %in% microbe_cols & to %in% vsc_targets) |
    (to %in% microbe_cols & from %in% vsc_targets)
  ) %>%
  mutate(
    microbe = ifelse(from %in% microbe_cols, from, to)
  ) %>%
  select(microbe) %>%
  distinct()

clr_long_all <- merged_data %>%
  filter(!is.na(diagnosis)) %>%
  pivot_longer(cols = all_of(microbe_cols),
               names_to = "microbe",
               values_to = "CLR_abundance") %>%
  filter(microbe %in% vsc_microbe_combined$microbe)

clr_summary_all <- clr_long_all %>%
  group_by(niche, diagnosis, microbe) %>%
  summarise(mean_abundance = mean(CLR_abundance, na.rm = TRUE), .groups = "drop")

clr_z_global <- clr_summary_all %>%
  group_by(microbe) %>%
  mutate(z_score = scale(mean_abundance)) %>%
  ungroup()

clr_z_global$niche <- factor(
  clr_z_global$niche,
  levels = c("interdental","subgingival","tongue")
)

clr_z_global$diagnosis <- factor(
  clr_z_global$diagnosis,
  levels = c("LPH","HPH","G","P")
)

# ------------------------------
# Heatmap plot
# ------------------------------

zscore_corr_plot <- ggplot(
  clr_z_global,
  aes(x = diagnosis, y = microbe, fill = z_score)
) +
  
  geom_tile(color = "white") +
  
  facet_wrap(
    ~ niche,
    scales = "free_y",
    nrow = 1,
    strip.position = "top",
    labeller = labeller(
      niche = c(
        interdental = "Interdental Plaque",
        subgingival = "Subgingival Plaque",
        tongue = "Tongue Dorsum"
      )
    )
  ) +
  
  scale_fill_distiller(
    palette = "RdBu",
    direction = -1,
    name = "Z-score"
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.placement = "outside"
  ) +
  
  labs(
    x = "Diagnosis Group",
    y = NULL
  )

# ------------------------------
# Save figure
# ------------------------------

ggsave(
  "figures/figure3/figure3_vsc_microbe_heatmap.svg",
  plot = zscore_corr_plot,
  width = 16,
  height = 10
)
