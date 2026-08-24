# R/R_snake/EA_42_get_umi_counts_in_umis_as_sublibs_mode.R

get_umi_counts_in_umis_as_sublibs_mode <- function(cfg){
  #-----------------------------------------------------------------------------
  # Read counts back into R
  #-----------------------------------------------------------------------------
  
  logger::log_info("Reading in read/UMI counts...")
  
  count_df_long <- process_bcwithqc_data(
    cfg = cfg,
    data_type = "umis"
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
    
    if (identical(combine_for_guide_stats, "sample")) {
      
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
  return(count_df_long)
}