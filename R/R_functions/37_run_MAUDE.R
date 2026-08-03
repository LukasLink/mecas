#-----------------------------------------------------------------------------
# Run MAUDE
#-----------------------------------------------------------------------------

run_MAUDE <- function(maude_counts_df,
                      cfg,
                      file_suffix = cfg$suffix$file_suffix,
                      run_maude_stage = TRUE,
                      run_plots_stage = FALSE) {
  
  maude_gene_stats_fpath <- file.path(
    cfg$paths$rds_output_folder,
    paste0("MAUDE_gene_stats", file_suffix)
  )
  
  maude_guide_stats_fpath <- file.path(
    cfg$paths$rds_output_folder,
    paste0("MAUDE_guide_stats", file_suffix)
  )
  
  unique_exp <- unique(maude_counts_df$exp)
  
  # Define binStats
  binStats <- get_binStats(bins = cfg$bins, experiments = unique_exp)
  
  sorted_bin_names <- cfg$bins %>%
    dplyr::filter(sorted_or_unsorted == "sorted") %>%
    dplyr::arrange(bin_fraction_min, bin_fraction_max) %>%
    dplyr::pull(bin_name)
  
  unsorted_bin_name <- cfg$bins %>%
    dplyr::filter(sorted_or_unsorted == "unsorted") %>%
    dplyr::pull(bin_name)
  
  if (isTRUE(run_maude_stage)) {
    
    maude_guide_stats <- findGuideHitsAllScreens(
      experiments = unique(maude_counts_df["exp"]),
      countDataFrame = maude_counts_df,
      binStats = binStats,
      sortBins = sorted_bin_names,
      unsortedBin = unsorted_bin_name,
      negativeControl = "isNontargeting"
    )
    
    logger::log_info("Saving MAUDE guide stats to: {maude_guide_stats_fpath}")
    saveRDS(maude_guide_stats, maude_guide_stats_fpath)
    
  } else if (isTRUE(run_plots_stage)) {
    
    .check_rds_exists(maude_guide_stats_fpath, "MAUDE guide stats")
    logger::log_info("Loading MAUDE guide stats from: {maude_guide_stats_fpath}")
    maude_guide_stats <- readRDS(maude_guide_stats_fpath)
    
  } else {
    stop_log("This should never trigger, check run_MAUDE.")
  }
  
  if (isTRUE(run_maude_stage)) {
    
    maude_guide_stats <- maude_guide_stats %>%
      dplyr::left_join(
        cfg$merged_sgRNA_df %>%
          dplyr::select(sgrna_id, entrez) %>%
          dplyr::distinct(sgrna_id, .keep_all = TRUE),
        by = c("sgRNA" = "sgrna_id")
      ) %>%
      dplyr::mutate(entrez = dplyr::coalesce(as.character(entrez), sgRNA))
    
    # Any entries from include_controls are manually turned into genes.
    include_controls_list <- cfg$controls$include_controls %||% character()
    
    if (length(include_controls_list) > 0) {
      for (control_gene in include_controls_list) {
        # Remove everything after the last _, so AAVS1_9 and AAVS1_13 are both AAVS1.
        control_gene <- sub("_[^_]*$", "", control_gene)
        maude_guide_stats$entrez[
          grepl(control_gene, maude_guide_stats$sgRNA)
        ] <- control_gene
      }
    }
    
    combine_for_gene_stats <- cfg$counting$combine_for_gene_stats %||% "none"
    
    if (!identical(combine_for_gene_stats, "none")) {
      if (!(combine_for_gene_stats %in% c("all", "sublib", "sample"))) {
        stop(
          "counting.combine_for_gene_stats must be one of: 'all', 'none', 'sublib', 'sample'",
          call. = FALSE
        )
      }
      
      if (identical(combine_for_gene_stats, "all")) {
        maude_guide_stats$exp <- "rep1"
      }
      
      if (identical(combine_for_gene_stats, "sublib")) {
        maude_guide_stats <- maude_guide_stats %>% 
          dplyr::mutate(exp = sample)
      }
      
      if (identical(combine_for_gene_stats, "sample")) {
        maude_guide_stats <- maude_guide_stats %>% 
          dplyr::mutate(exp = sublib)
      }
    }
    
    maude_gene_stats <- getElementwiseStats(
      experiments = unique(maude_guide_stats["exp"]),
      normNBSummaries = maude_guide_stats,
      negativeControl = "isNontargeting",
      elementIDs = "entrez"
    )
    
    maude_gene_stats <- maude_gene_stats %>%
      dplyr::filter(numGuides >= cfg$filtering$min_guides_per_gene)
    
    logger::log_info("Saving MAUDE gene stats to: {maude_gene_stats_fpath}")
    saveRDS(maude_gene_stats, maude_gene_stats_fpath)
    
  } else if (isTRUE(run_plots_stage)) {
    
    .check_rds_exists(maude_gene_stats_fpath, "MAUDE gene stats")
    logger::log_info("Loading MAUDE gene stats from: {maude_gene_stats_fpath}")
    maude_gene_stats <- readRDS(maude_gene_stats_fpath)
    
  } else {
    stop_log("This should never trigger, check start_with processing.")
  }
  
  maude_guide_stats <- readRDS(maude_guide_stats_fpath)
  maude_gene_stats <- readRDS(maude_gene_stats_fpath)
  
  invisible(list(
    maude_guide_stats = maude_guide_stats,
    maude_gene_stats = maude_gene_stats,
    maude_bins = binStats,
    maude_guide_stats_fpath = maude_guide_stats_fpath,
    maude_gene_stats_fpath = maude_gene_stats_fpath
  ))
}

.check_rds_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop_log(
      label, " file does not exist:\n  ",
      path,
      "\n\nSet `first_time: true` to generate it, or check your output/rds folder and suffix settings."
    )
  }
  
  invisible(TRUE)
}