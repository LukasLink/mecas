# R/R_functions/CA_32_make_reads_per_umi_coverage_df.R

make_reads_per_umi_coverage_df <- function(count_df_long, cfg) {
  
  log_info("Coverage start | RSS: {sprintf('%.2f GB', memory_rss_gb())}")
  
  umi_per_guide <- count_df_long %>%
    dplyr::distinct(
      bin_name,
      sample,
      sgRNA,
      sublib
    ) %>%
    dplyr::count(
      bin_name,
      sample,
      sgRNA,
      name = "n_umis"
    )
  log_info("After umi_per_guide | RSS: {sprintf('%.2f GB', memory_rss_gb())}")
  
  total_guides <- dplyr::n_distinct(cfg$library$sgrna_id)
  
  guide_coverage <- umi_per_guide %>%
    dplyr::group_by(bin_name, sample) %>%
    dplyr::summarise(
      guides_observed = dplyr::n_distinct(sgRNA),
      coverage_percent = round(guides_observed / total_guides * 100, digits = 3),
      total_umi_sub_libraries = sum(n_umis),
      mean_umis_per_guide = mean(n_umis, na.rm = TRUE),
      median_umis_per_guide = median(n_umis, na.rm = TRUE),
      sd_umis_per_guide = sd(n_umis, na.rm = TRUE),
      min_umis_per_guide = min(n_umis, na.rm = TRUE),
      max_umis_per_guide = max(n_umis, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  log_info("After guide_coverage | RSS: {sprintf('%.2f GB', memory_rss_gb())}")
  
  reads_per_umi <- count_df_long %>%
    dplyr::group_by(bin_name, sample) %>%
    dplyr::summarise(
      mean_reads_per_umi = mean(count, na.rm = TRUE),
      median_reads_per_umi = median(count, na.rm = TRUE),
      sd_reads_per_umi = sd(count, na.rm = TRUE),
      min_reads_per_umi = min(count, na.rm = TRUE),
      max_reads_per_umi = max(count, na.rm = TRUE),
      total_reads = sum(count, na.rm = TRUE),
      .groups = "drop"
    )
  
  log_info("After reads_per_umi | RSS: {sprintf('%.2f GB', memory_rss_gb())}")
  
  guide_coverage %>%
    dplyr::left_join(
      reads_per_umi,
      by = c("bin_name", "sample")
    )
}