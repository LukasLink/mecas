# R/R_functions/CA_20_get_read_count_df_long.R

get_read_count_df_long <- function(cfg){
  #-----------------------------------------------------------------------------
  # Read counts back into R
  #-----------------------------------------------------------------------------
  
  
  count_df_long <- process_bcwithqc_data(
    cfg = cfg,
    data_type = "reads"
  )

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
    pre_combining_fpath <- file.path(cfg$paths$rds_output_folder, "raw_reads_count_df_long_pre_combining.rds")
    
    if (identical(combine_for_guide_stats, "sample")) {
      log_info("Combining Samples for MAUDE guide stat calculation. Saving pre-combining raw read counts to {pre_combining_fpath}")
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
      log_info("Combining Sublibraries for MAUDE guide stat calculation. Saving pre-combining raw read counts to {pre_combining_fpath}")
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
  count_df_tsv_gz <- file.path(
    cfg$paths$results_output_folder,
    "raw_reads_count_df.tsv.gz"
  )
  
  # R-specific compact version
  saveRDS(
    count_df_long,
    file = cfg$paths$reads_count_df_fpath,
    compress = "gzip"
  )
  
  # Compressed, portable TSV version
  readr::write_tsv(
    count_df_long,
    file = count_df_tsv_gz
  )
}