# R/R_functions/CA_33_apply_optional_count_df_operations.R


apply_optional_count_df_operations <- function(count_df_long, cfg){
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
      
      log_info(
        "Before combining samples | R memory: {sprintf('%.2f GB', memory_used_gb())} | RSS: {sprintf('%.2f GB', memory_rss_gb())}"
      )
      # Using data.table here because this operation balooned R memory requirements
      data.table::setDT(count_df_long)
      
      count_df_long <- count_df_long[
        ,
        .(
          count = sum(count, na.rm = TRUE),
          group_category = data.table::first(group_category)
        ),
        by = .(sgRNA, sublib, bin_name)
      ]
      
      count_df_long[
        ,
        `:=`(
          sample = "combined_samples",
          exp = paste(sublib, "combined_samples", sep = "_")
        )
      ]
    }
    log_info(
      "After combining samples | R memory: {sprintf('%.2f GB', memory_used_gb())} | RSS: {sprintf('%.2f GB', memory_rss_gb())}"
    )
    
    if (identical(combine_for_guide_stats, "sublib")) {
      log_info("Combining Sublibraries for MAUDE guide stat calculation. Saving pre-combining counts to {pre_combining_fpath}")
      saveRDS(count_df_long, pre_combining_fpath)
      
      # Using data.table here to reduce memory requirements for large count dataframes.
      data.table::setDT(count_df_long)
      
      count_df_long <- count_df_long[
        ,
        .(
          count = sum(count, na.rm = TRUE),
          group_category = data.table::first(group_category)
        ),
        by = .(sgRNA, sample, bin_name)
      ]
      
      count_df_long[
        ,
        `:=`(
          sublib = "combined_sublibraries",
          exp = paste("combined_sublibraries", sample, sep = "_")
        )
      ]
    }
    log_info("Finished Combining.")
  }
  
  return(count_df_long)
}