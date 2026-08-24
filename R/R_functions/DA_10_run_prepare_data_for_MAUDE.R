# R/R_functions/DA_10_run_prepare_data_for_MAUDE.R

#-----------------------------------------------------------------------------
# prepare data for MAUDE
#-----------------------------------------------------------------------------

run_prepare_data_for_MAUDE <- function(count_df_long,
                                       cfg) {
  logger::log_info("Preparing Data for MAUDE...")
  
  pseudocount_already_added <- FALSE
  check_all_bins_present(count_df_long = count_df_long, cfg = cfg)
  
  binStats_fpath <- file.path(cfg$paths$rds_output_folder, "MAUDE_binStats.rds")
  saveRDS(cfg$bins, binStats_fpath)
  logger::log_info("Bins validated and saved to {binStats_fpath}")
  #-----------------------------------------------------------------------------
  # Optional: subsample controls
  #-----------------------------------------------------------------------------
  
  if (isTRUE(cfg$controls$subsample_controls)) {
    count_df_long <- subsample_controls_func(
      count_df_long = count_df_long,
      cfg = cfg
    )
  }

  #-----------------------------------------------------------------------------
  # Optional: normalize counts
  #-----------------------------------------------------------------------------
  
  norm_method <- cfg$normalization$norm_method %||% ""
  
  if (!(norm_method %in% c("", "control_median"))) {
    stop(
      "Error: `normalization.norm_method` must be one of '', 'control_median'.",
      call. = FALSE
    )
  }
  log_info(
    "before normalize_count_df_long | R memory: {sprintf('%.2f GB', memory_used_gb())} | RSS: {sprintf('%.2f GB', memory_rss_gb())}"
  )
  if (identical(norm_method, "control_median")) {
    list2env(
      normalize_count_df_long(
        count_df_long = count_df_long,
        umis_as_sublibs = cfg$counting$umis_as_sublibs,
        norm_method = norm_method,
        return_info = TRUE
      ),
      envir = environment()
    )
  }
  log_info(
    "After normalize_count_df_long | R memory: {sprintf('%.2f GB', memory_used_gb())} | RSS: {sprintf('%.2f GB', memory_rss_gb())}"
  )
  #-----------------------------------------------------------------------------
  # Convert long count table to MAUDE wide format
  #-----------------------------------------------------------------------------
  
  maude_counts_df <- count_df_long_to_wide(
    count_df_long = count_df_long,
    bins = cfg$bins,
    print = FALSE,
    drop_0s = cfg$filtering$drop_0s,
    pseudocount_already_added = pseudocount_already_added,
    recover_input = cfg$normalization$recover_input,
    strict_mode = cfg$filtering$strict_mode,
    umis_as_sublibs = cfg$counting$umis_as_sublibs
  )
  rm(count_df_long)
  log_info(
    "After count_df_long_to_wide | R memory: {sprintf('%.2f GB', memory_used_gb())} | RSS: {sprintf('%.2f GB', memory_rss_gb())}"
  )
  #-----------------------------------------------------------------------------
  # Replicate handling
  #-----------------------------------------------------------------------------
  
  method <- cfg$replicates$method %||% ""
  
  if (!(method %in% c("", "rep", "rep_sample", "rep_sublib"))) {
    stop(
      "Error: `replicates.method` must be one of '', 'rep', 'rep_sample', or 'rep_sublib'. ",
      "Method 'sum' is deprecated; use `counting.combine_for_guide_stats` instead.",
      call. = FALSE
    )
  }
  
  if (identical(method, "")) {
    maude_counts_df <- maude_counts_df %>%
      dplyr::mutate(exp = "rep1")
  }
  
  if (identical(method, "rep")) {
    # Keep existing exp column.
  }
  
  if (identical(method, "rep_sample")) {
    maude_counts_df <- maude_counts_df %>%
      dplyr::mutate(exp = sample)
  }
  
  if (identical(method, "rep_sublib")) {
    maude_counts_df <- maude_counts_df %>%
      dplyr::mutate(exp = sublib)
  }
  

  
  #-----------------------------------------------------------------------------
  # Optional: control handling
  #-----------------------------------------------------------------------------
  
  include_controls_list <- cfg$controls$include_controls %||% character()
  use_only_these_controls_list <- cfg$controls$use_only_these_controls %||% character()
  
  if (length(include_controls_list) > 0) {
    for (control_gene in include_controls_list) {
      maude_counts_df$isNontargeting[
        grepl(control_gene, maude_counts_df$sgRNA)
      ] <- FALSE
    }
  }
  
  if (length(use_only_these_controls_list) > 0) {
    maude_counts_df$isNontargeting[
      !(maude_counts_df$sgRNA %in% use_only_these_controls_list)
    ] <- FALSE
  }
  
  logger::log_info("Finished preparing Data for MAUDE.")
  
  saveRDS(maude_counts_df, cfg$paths$maude_counts_df_fpath)
  logger::log_info("Saving MAUDE ready dataframe to: {cfg$paths$maude_counts_df_fpath}")
  
  return(maude_counts_df)
}