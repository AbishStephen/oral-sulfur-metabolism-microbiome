# ========================================
# Figure 4
# Microbe–cytokine correlations
# ========================================

library(tidyverse)
library(ggplot2)
library(psych)
library(compositions)

sample_info <- read.csv("data/sampleinfo.csv")
clinical_data <- read.csv("data/clinical_data.csv")
gcf_cytokines <- read.csv("data/gcf_cytokines.csv")
microbial_data <- read.csv("data/microbial_data.csv")

sample_info <- sample_info[-c(4)]

microbial_with_info <- microbial_data %>%
  left_join(sample_info, by = "sample")

merged_data <- microbial_with_info %>%
  left_join(clinical_data, by = "patno") %>%
  left_join(gcf_cytokines, by = "patno")

merged_data <- merged_data %>%
  mutate(
    h2s_log = log1p(h2s),
    ch3sh_log = log1p(ch3sh),
    gcf_h2s_um_log = log1p(gcf_h2s_um),
    gcf_ch3sh_um_log = log1p(gcf_ch3sh_um)
  )

cytokine_cols <- names(gcf_cytokines)[!(names(gcf_cytokines) %in% c("patno"))]
microbe_cols <- setdiff(names(microbial_data), "sample")

merged_data[cytokine_cols] <- scale(merged_data[cytokine_cols])
merged_data[microbe_cols] <- merged_data[microbe_cols] + 1
merged_data[microbe_cols] <- clr(merged_data[microbe_cols])

filtered_data <- merged_data %>%
  filter(diagnosis %in% c("LPH","HPH","G"))

niches <- c("interdental","subgingival","tongue")

significant_correlations_cytokine_microbe <- list()

for (niche in niches) {

  niche_data <- filtered_data %>%
    filter(niche == !!niche)

  selected_data <- niche_data %>%
    select(
      h2s_log,
      ch3sh_log,
      gcf_h2s_um_log,
      gcf_ch3sh_um_log,
      all_of(cytokine_cols),
      all_of(microbe_cols)
    )

  if (nrow(selected_data) > 10) {

    cor_matrix <- corr.test(selected_data, method = "spearman")

    cor_coeffs <- cor_matrix$r
    p_values <- cor_matrix$p

    p_adj <- p.adjust(p_values, method = "fdr")

    cor_coeffs_filtered <- cor_coeffs
    cor_coeffs_filtered[p_adj > 0.05 | abs(cor_coeffs_filtered) < 0.4] <- NA

    sig_pairs <- which(!is.na(cor_coeffs_filtered), arr.ind = TRUE)

    sig_correlations <- data.frame(
      from = rownames(cor_coeffs_filtered)[sig_pairs[,1]],
      to = colnames(cor_coeffs_filtered)[sig_pairs[,2]],
      weight = cor_coeffs_filtered[sig_pairs]
    )

    sig_correlations <- sig_correlations %>%
      filter(
        (from %in% microbe_cols & to %in% cytokine_cols) |
        (from %in% cytokine_cols & to %in% microbe_cols)
      ) %>%
      mutate(niche = niche)

    significant_correlations_cytokine_microbe[[niche]] <- sig_correlations
  }
}

dot_data <- significant_correlations_cytokine_microbe[["interdental"]]

dot_data_filtered <- dot_data %>%
  mutate(
    microbe = ifelse(from %in% microbe_cols, from, to),
    cytokine = ifelse(from %in% cytokine_cols, from, to),
    corr_sign = ifelse(weight > 0, "Positive", "Negative"),
    abs_weight = abs(weight)
  )

interdental_dotheatmap <- ggplot(dot_data_filtered,
                                 aes(x = cytokine, y = microbe)) +
  geom_point(aes(size = abs_weight, color = corr_sign)) +
  theme_minimal()

ggsave(
  "figures/figure4/interdental_dotheatmap.svg",
  interdental_dotheatmap,
  height = 15,
  width = 10
)

#subgingival
dot_data <- significant_correlations_cytokine_microbe[["subgingival"]]

# Filter microbes and cytokines 
dot_data_filtered <- dot_data %>%
  filter((from %in% microbe_cols & to %in% cytokine_cols) |
           (from %in% cytokine_cols & to %in% microbe_cols)) %>%
  mutate(
    microbe = ifelse(from %in% microbe_cols, from, to),
    cytokine = ifelse(from %in% cytokine_cols, from, to),
    corr_sign = ifelse(weight > 0, "Positive", "Negative"),
    abs_weight = abs(weight)
  )

# Create the dot heatmap
subgingival_dotheatmap <- ggplot(dot_data_filtered, aes(x = cytokine, y = microbe)) +
  geom_point(aes(size = abs_weight, color = corr_sign)) +
  scale_color_manual(values = c("Positive" = "firebrick", "Negative" = "steelblue")) +
  scale_size(range = c(2, 6)) +
  theme_minimal() +
  labs(title = "Subgingival Plaque",
       x = "Cytokine",
       y = "Microbial Taxa",
       size = "|Correlation Coefficient|",
       color = "Correlation Sign") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("subgingival_dotheatmap.svg", plot = subgingival_dotheatmap, height = 15, width = 10)

#tongue
dot_data <- significant_correlations_cytokine_microbe[["tongue"]]

# Filter microbes and cytokines 
dot_data_filtered <- dot_data %>%
  filter((from %in% microbe_cols & to %in% cytokine_cols) |
           (from %in% cytokine_cols & to %in% microbe_cols)) %>%
  mutate(
    microbe = ifelse(from %in% microbe_cols, from, to),
    cytokine = ifelse(from %in% cytokine_cols, from, to),
    corr_sign = ifelse(weight > 0, "Positive", "Negative"),
    abs_weight = abs(weight)
  )

# Create the dot heatmap
tongue_dotheatmap <- ggplot(dot_data_filtered, aes(x = cytokine, y = microbe)) +
  geom_point(aes(size = abs_weight, color = corr_sign)) +
  scale_color_manual(values = c("Positive" = "firebrick", "Negative" = "steelblue")) +
  scale_size(range = c(2, 6)) +
  theme_minimal() +
  labs(title = "Tongue",
       x = "Cytokine",
       y = "Microbial Taxa",
       size = "|Correlation Coefficient|",
       color = "Correlation Sign") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("tongue_dotheatmap.svg", plot = tongue_dotheatmap, height = 15, width = 10)
