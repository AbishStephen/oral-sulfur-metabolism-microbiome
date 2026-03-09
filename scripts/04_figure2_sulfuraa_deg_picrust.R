# ========================================
# Figure 2
# Sulfur metabolism and cysteine degradation
# ========================================

library(tidyverse)
library(reshape2)
library(ggplot2)
library(psych)
library(patchwork)
library(compositions)

# ------------------------------
# Load data
# ------------------------------

sample_info <- read.csv("data/sampleinfo.csv")
clinical_data <- read.csv("data/clinical_data.csv")
gcf_cytokines <- read.csv("data/gcf_cytokines.csv")
microbial_data <- read.csv("data/microbial_data.csv")

sample_info <- sample_info[-c(4)]

# ------------------------------
# Merge datasets
# ------------------------------

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

# ------------------------------
# PICRUSt pathway analysis
# ------------------------------

picrust_pathways <- read.table(
  "data/picrust/path_abun_unstrat.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

picrust_pathways <- as.data.frame(t(picrust_pathways))
picrust_pathways$sample <- rownames(picrust_pathways)

picrust_merged <- picrust_pathways %>%
  left_join(sample_info, by = "sample") %>%
  left_join(clinical_data, by = "patno")

meta_cols <- c(
  "sample","patno","niche","diagnosis","h2s","ch3sh",
  "gcf_h2s_um","gcf_ch3sh_um"
)

pathway_cols <- setdiff(names(picrust_merged), meta_cols)

picrust_long <- picrust_merged %>%
  pivot_longer(
    cols = all_of(pathway_cols),
    names_to = "pathway",
    values_to = "abundance"
  )

# ------------------------------
# Sulfur pathway filtering
# ------------------------------

sulfur_pathways <- picrust_long %>%
  filter(str_detect(tolower(pathway), "cysteine|methionine|sulfur|sulphur|sulfide"))

sulfur_summary <- sulfur_pathways %>%
  group_by(niche, diagnosis, pathway) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE)) %>%
  ungroup()

sulfur_summary$diagnosis <- factor(
  sulfur_summary$diagnosis,
  levels = c("LPH","HPH","G","P"),
  labels = c("Low Plaque Health","High Plaque Health","Gingivitis","Periodontitis")
)

plot_list <- list()

for (niche_name in unique(sulfur_summary$niche)) {

  niche_data <- sulfur_summary %>% filter(niche == niche_name)

  p <- ggplot(
    niche_data,
    aes(
      x = mean_abundance + 1,
      y = pathway,
      color = diagnosis
    )
  ) +
    geom_point(size = 3) +
    scale_x_continuous(trans = "log10") +
    theme_minimal()

  plot_list[[niche_name]] <- p
}

final_panel <- wrap_plots(plot_list, ncol = 1)

ggsave(
  "figures/figure2/sulfur_metabolism_panel.svg",
  final_panel,
  width = 10,
  height = 18
)
