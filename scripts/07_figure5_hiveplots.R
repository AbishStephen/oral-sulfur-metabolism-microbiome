# ========================================
# Figure 5
# Microbe–cytokine–VSC networks
# ========================================

library(tidyverse)
library(igraph)
library(ggraph)
library(psych)
library(patchwork)
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

merged_data$diagnosis <- recode(
  merged_data$diagnosis,
  "LPH" = "H",
  "HPH" = "H"
)

cytokine_cols <- names(gcf_cytokines)[!(names(gcf_cytokines) %in% c("patno"))]
microbe_cols <- setdiff(names(microbial_data), "sample")

merged_data[cytokine_cols] <- scale(merged_data[cytokine_cols])
merged_data[microbe_cols] <- merged_data[microbe_cols] + 1
merged_data[microbe_cols] <- clr(merged_data[microbe_cols])

niches <- c("interdental","subgingival","tongue")

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
        from = rownames(cor_coeffs_filtered)[sig_pairs[,1]],
        to = colnames(cor_coeffs_filtered)[sig_pairs[,2]],
        weight = cor_coeffs_filtered[sig_pairs]
      )

      sig_correlations <- sig_correlations %>%
        filter(from != to) %>%
        mutate(niche = niche, diagnosis = diagnosis)

      significant_correlations[[paste(niche, diagnosis, sep = "_")]] <- sig_correlations
    }
  }
}

