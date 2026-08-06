# R/R_functions_snake/11_get_count_df_long_and_make_coverage_file.R

get_count_df_long_and_make_coverage_file <- function(cfg){
  #-----------------------------------------------------------------------------
  # Read counts back into R
  #-----------------------------------------------------------------------------
  
  logger::log_info("Reading in read/UMI counts...")
  
  count_df_long <- process_bcwithqc_data(
    cfg = cfg,
    data_type = cfg$counting$data_type
  )
  
  logger::log_info("Finished reading in read/UMI counts.")
  logger::log_info("sgRNAs aligned to the wrong sublibrary were excluded from the analysis.")
  
  #-----------------------------------------------------------------------------
  # Optional: use_only_these_controls / include_controls
  #-----------------------------------------------------------------------------
  
  use_only_these_controls_list <- cfg$controls$use_only_these_controls %||% character()
  include_controls_list <- cfg$controls$include_controls %||% character()
  
  if (length(use_only_these_controls_list) > 0) {
    excluded_controls <- count_df_long %>%
      dplyr::filter(
        group_category != "targeting",
        !sgRNA %in% use_only_these_controls_list
      ) %>%
      dplyr::distinct(sgRNA) %>%
      dplyr::pull(sgRNA)
    
    include_controls_list <- c(include_controls_list, excluded_controls)
  }
  
  if (length(include_controls_list) > 0) {
    for (control_gene in include_controls_list) {
      count_df_long$group_category[
        grepl(control_gene, count_df_long$sgRNA)
      ] <- "targeting"
    }
  }
  
  #-----------------------------------------------------------------------------
  # Optional: combine samples or sublibraries of the same bin_name
  #-----------------------------------------------------------------------------
  
  combine_for_guide_stats <- cfg$replicates$combine_for_guide_stats %||% ""
  
  if (!identical(combine_for_guide_stats, "")) {
    pre_combining_fpath <- file.path(cfg$paths$rds_output_folder, "count_df_long_pre_combining.rds")
    
    if (identical(combine_for_guide_stats, "sample")) {
      log_info("Combining Samples for MAUDE guide stat calculation. Saving pre-combining counts to {pre_combining_fpath}")
      saveRDS(count_df_long, pre_combining_fpath)
      
      count_df_long <- count_df_long %>%
        dplyr::group_by(sgRNA, sublib, bin_name) %>%
        dplyr::summarise(
          count = sum(count, na.rm = TRUE),
          group_category = dplyr::first(group_category),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          sample = "combined_samples",
          exp = paste(sublib, sample, sep = "_")
        )
    }
    
    if (identical(combine_for_guide_stats, "sublib")) {
      log_info("Combining Sublibraries for MAUDE guide stat calculation. Saving pre-combining counts to {pre_combining_fpath}")
      saveRDS(count_df_long, pre_combining_fpath)
      
      count_df_long <- count_df_long %>%
        dplyr::group_by(sgRNA, sample, bin_name) %>%
        dplyr::summarise(
          count = sum(count, na.rm = TRUE),
          group_category = dplyr::first(group_category),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          sublib = "combined_sublibraries",
          exp = paste(sublib, sample, sep = "_")
        )
    }
  }
  
  #-----------------------------------------------------------------------------
  # Make coverage file
  #-----------------------------------------------------------------------------
  
  logger::log_info("Creating coverage file for read/UMI counts...")
  
  mapping_results_df <- parse_coverage_file(cfg$paths$log_files$count)

  overall_targeting <- count_df_long %>%
    dplyr::summarise(
      total_counts = sum(count, na.rm = TRUE),
      targeting_counts = sum(count[group_category == "targeting"], na.rm = TRUE),
      targeting_perc = 100 * targeting_counts / total_counts
    ) %>%
    dplyr::mutate(targeting_perc = sprintf("%.2f%%", targeting_perc))
  
  targeting_by_group <- count_df_long %>%
    dplyr::group_by(bin_name, sublib, sample) %>%
    dplyr::summarise(
      total_counts = sum(count, na.rm = TRUE),
      targeting_counts = sum(count[group_category == "targeting"], na.rm = TRUE),
      targeting_perc = 100 * targeting_counts / total_counts,
      .groups = "drop"
    ) %>%
    dplyr::mutate(targeting_perc = sprintf("%.2f%%", targeting_perc))
  
  # overall_targeting_merged <- cfg$merged_sgRNA_df %>%
  #   dplyr::summarise(
  #     total_counts = sum(count, na.rm = TRUE),
  #     targeting_counts = sum(count[!is.na(entrez)], na.rm = TRUE),
  #     targeting_perc = 100 * targeting_counts / total_counts
  #   ) %>%
  #   dplyr::mutate(targeting_perc = sprintf("%.2f%%", targeting_perc))
  
  overall_targeting_merged <- data.frame(
    note = "overall_targeting_reference temporarily disabled in Snakemake mode"
  )
  
  writexl::write_xlsx(
    list(
      mapping_results = mapping_results_df,
      overall_targeting = overall_targeting,
      targeting_by_group = targeting_by_group,
      overall_targeting_reference = overall_targeting_merged
    ),
    file.path(cfg$paths$results_output_folder, "mapping_results.xlsx")
  )
  
  logger::log_info("Finished creating coverage file for read/UMI counts.")
  logger::log_info(
    paste(
      "The coverage file is:",
      file.path(cfg$paths$results_output_folder, "mapping_results.xlsx")
    )
  )
  
  return(count_df_long)
}