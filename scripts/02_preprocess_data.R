library(compositions)

sample_info <- sample_info[-c(4)]

microbial_with_info <- microbial_data %>%
  left_join(sample_info, by = "sample")

merged_data <- microbial_with_info %>%
  left_join(clinical_data, by = "patient_no") %>%
  left_join(gcf_cytokines, by = "patient_no")

merged_data <- merged_data %>%
  mutate(
    h2s_log = log1p(h2s),
    ch3sh_log = log1p(ch3sh),
    gcf_h2s_um_log = log1p(gcf_h2s_um),
    gcf_ch3sh_um_log = log1p(gcf_ch3sh_um)
  )

cytokine_cols <- names(gcf_cytokines)[!(names(gcf_cytokines) %in% c("patient_no"))]

microbe_cols <- setdiff(names(microbial_data), "sample")

merged_data[cytokine_cols] <- scale(merged_data[cytokine_cols])

merged_data[microbe_cols] <- merged_data[microbe_cols] + 1
merged_data[microbe_cols] <- clr(merged_data[microbe_cols])
