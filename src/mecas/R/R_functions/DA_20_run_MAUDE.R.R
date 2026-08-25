# R/R_functions/DA_20_run_MAUDE.R

#-----------------------------------------------------------------------------
# Run MAUDE
#-----------------------------------------------------------------------------

run_MAUDE <- function(maude_counts_df,
                      cfg,
                      file_suffix = cfg$suffix$file_suffix,
                      run_maude_stage = TRUE,
                      run_plots_stage = FALSE) {
  
  
  # TESTING ONLY REMOVE THIS IN FINAL VERSION
  #############################################################################
  # set.seed(123)
  # test_pre_guide <- maude_counts_df
  # 
  # # Add random background variation of 1-3 counts to every bin.
  # bin_cols <- c("I", "L", "U")
  # 
  # for (bin_col in bin_cols) {
  #   test_pre_guide[[bin_col]] <- test_pre_guide[[bin_col]] +
  #     sample(1:3, nrow(test_pre_guide), replace = TRUE)
  # }
  # 
  # # Add a strong positive effect to U for selected genes.
  # upper_genes <- c("ABCA1", "ACSS3", "AAMP")
  # 
  # upper_rows <- grepl(
  #   paste0("^(", paste(upper_genes, collapse = "|"), ")_[0-9]+$"),
  #   test_pre_guide$sgRNA
  # )
  # 
  # test_pre_guide$U[upper_rows] <-
  #   test_pre_guide$U[upper_rows] +
  #   sample(10:20, sum(upper_rows), replace = TRUE)
  # 
  # # Add a strong negative effect to L for selected genes.
  # lower_genes <- c("ZP3", "ZNF710", "ZCWPW2")
  # 
  # lower_rows <- grepl(
  #   paste0("^(", paste(lower_genes, collapse = "|"), ")_[0-9]+$"),
  #   test_pre_guide$sgRNA
  # )
  # 
  # test_pre_guide$L[lower_rows] <-
  #   test_pre_guide$L[lower_rows] +
  #   sample(10:20, sum(lower_rows), replace = TRUE)
  # 
  # maude_counts_df <- test_pre_guide
  # rm(test_pre_guide)
  # 
  # log_warn("TESTING ENABLED FALSE DATA WILL BE USED!!!!!!!!")
  #############################################################################
  
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
    if (isTRUE(cfg$debug$MAUDE)){log_info("Before Guide stats | R memory: {sprintf('%.2f GB', memory_used_gb())} | RSS: {sprintf('%.2f GB', memory_rss_gb())}")
    maude_guide_stats <- findGuideHitsAllScreens(
      experiments = unique(maude_counts_df["exp"]),
      countDataFrame = maude_counts_df,
      binStats = binStats,
      sortBins = sorted_bin_names,
      unsortedBin = unsorted_bin_name,
      negativeControl = "isNontargeting"
    )
    if (isTRUE(cfg$debug$MAUDE)){log_info("After Guide stats | R memory: {sprintf('%.2f GB', memory_used_gb())} | RSS: {sprintf('%.2f GB', memory_rss_gb())}")
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
    # --------------------------------------------------------------------------
    # Check for degenerate guide statistics before gene-level MAUDE analysis.
    # --------------------------------------------------------------------------
    
    control_z <- maude_guide_stats$Z[
      maude_guide_stats$isNontargeting %in% TRUE
    ]
    
    control_z <- control_z[is.finite(control_z)]
    
    if (
      length(control_z) < 2 ||
      isTRUE(all(control_z == control_z[1])) ||
      isTRUE(sd(control_z) == 0)
    ) {
      debug_pre_guide_fpath <- file.path(
        cfg$paths$rds_output_folder,
        "DEBUG_maude_df_before_guide_stats.rds"
      )
      
      debug_post_guide_fpath <- file.path(
        cfg$paths$rds_output_folder,
        "DEBUG_maude_df_after_guide_stats.rds"
      )
      
      saveRDS(maude_counts_df, debug_pre_guide_fpath)
      saveRDS(maude_guide_stats, debug_post_guide_fpath)
      
      stop_log(
        "MAUDE gene-level statistics cannot be calculated because the ",
        "non-targeting guide Z scores have no variance.\n\n",
        "This usually means that the screen contains no measurable differences ",
        "between the sorted bins. This can occur, for example, when identical ",
        "or effectively identical count data are provided for multiple bins.\n\n",
        "Number of finite non-targeting guide Z scores: ", length(control_z), "\n",
        "Non-targeting guide Z-score SD: ",
        ifelse(length(control_z) >= 2, sd(control_z), NA_real_), "\n\n",
        "Debugging dataframes were saved to:\n",
        "  Before guide statistics: ", debug_pre_guide_fpath, "\n",
        "  After guide statistics:  ", debug_post_guide_fpath
      )
    }
    # --------------------------------------------------------------------------
    if (isTRUE(cfg$debug$MAUDE)){log_info("Before Gene stats | R memory: {sprintf('%.2f GB', memory_used_gb())} | RSS: {sprintf('%.2f GB', memory_rss_gb())}")
    maude_gene_stats <- getElementwiseStats(
      experiments = unique(maude_guide_stats["exp"]),
      normNBSummaries = maude_guide_stats,
      negativeControl = "isNontargeting",
      elementIDs = "entrez"
    )
    if (isTRUE(cfg$debug$MAUDE)){log_info("After Gene stats | R memory: {sprintf('%.2f GB', memory_used_gb())} | RSS: {sprintf('%.2f GB', memory_rss_gb())}")
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