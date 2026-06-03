#-----------------------------------------------------------------------------
# Run MAUDE
#-----------------------------------------------------------------------------
run_MAUDE <- function(maude_counts_df, opt, file_suffix,
                      run_maude_stage,
                      run_plots_stage){
  
  maude_gene_stats_fpath <- file.path(rds_output_folder,
                                      paste0("MAUDE_gene_stats", file_suffix))
  maude_guide_stats_fpath <- file.path(rds_output_folder,
                                       paste0("MAUDE_guide_stats", file_suffix))
  
  unique_exp <- unique(maude_counts_df$exp)
  
  # Define bin stats
  lower_bin_end = upper_lower_percentage
  upper_bin_start = 1 - upper_lower_percentage
  
  maude_bins <- tibble(Bin = rep(c('upper', 'lower'), length(unique_exp)),  # Repeat 'upper' and 'lower' for each exp
                       exp = rep(unique_exp, each = 2),  # Repeat each exp value twice for 'upper' and 'lower'
                       binStartQ = ifelse(rep(c('upper', 'lower'), length(unique_exp)) == 'lower', 0.001, upper_bin_start),
                       binEndQ = ifelse(rep(c('upper', 'lower'), length(unique_exp)) == 'lower', lower_bin_end, 0.999),
                       fraction = binEndQ - binStartQ,
                       binStartZ = qnorm(binStartQ),
                       binEndZ = qnorm(binEndQ)) %>%
    select(Bin, binStartQ, binEndQ, fraction, binStartZ, binEndZ, exp) %>%
    as.data.frame()
  
  if (isTRUE(run_maude_stage)) {


    maude_guide_stats <- findGuideHitsAllScreens(
      experiments = unique(maude_counts_df["exp"]),
      countDataFrame = maude_counts_df,
      binStats = maude_bins,
      sortBins = c("lower", "upper"),
      unsortedBin = "input",
      negativeControl = "isNontargeting"
    )
    
    logger::log_info("Saving MAUDE guide stats to: {maude_guide_stats_fpath}")
    saveRDS(maude_guide_stats, maude_guide_stats_fpath)
    
  } else if (isTRUE(run_plots_stage)) {
    .check_rds_exists(maude_guide_stats_fpath, "MAUDE guide stats")
    logger::log_info("Loading MAUDE guide stats from: {maude_guide_stats_fpath}")
    maude_guide_stats <- readRDS(maude_guide_stats_fpath)
  } else {
    stop_log("This should never trigger, check start_with processing.")
  }
  
  
  if (isTRUE(run_maude_stage)){
    
    maude_guide_stats <- maude_guide_stats %>%
      left_join(
        merged_sgRNA_df %>%
          select(sgrna_id, entrez) %>%
          distinct(sgrna_id, .keep_all = TRUE),
        by = c("sgRNA" = "sgrna_id")
      ) %>%
      mutate(entrez = coalesce(as.character(entrez), sgRNA))
    
    # any entries from include_controls_list are manually turned into genes
    if (length(include_controls_list) > 0) {
      for (control_gene in include_controls_list) {
        # here we remove everything after the last _, so stuff like AAVS1_9 and
        # AAVS1_13 are both treated as AAVS1
        control_gene <- sub("_[^_]*$", "", control_gene)
        maude_guide_stats$entrez[grepl(control_gene, maude_guide_stats$sgRNA)] <- control_gene
      }
    }
    
    if (combine_for_gene_stats != "none"){
      if (!(combine_for_gene_stats %in% c("all","sublib","sample"))){
        stop("combine_for_gene_stats must be one of: 'all','none','sublib','sample'")
      }
      if (combine_for_gene_stats == "all"){
        maude_guide_stats$exp <- "rep1"
      }
      if (combine_for_gene_stats == "sublib"){
        maude_guide_stats <- maude_guide_stats %>% 
          mutate(exp = sample)
      }
      if (combine_for_gene_stats == "sample"){
        maude_guide_stats <- maude_guide_stats %>% 
          mutate(exp = sublib)
      }    
    }
    
    ## calculate gene-level summarized scores
    maude_gene_stats <- getElementwiseStats(
      experiments = unique(maude_guide_stats['exp']),
      normNBSummaries = maude_guide_stats,
      negativeControl = 'isNontargeting',
      elementIDs = 'entrez'
    )
    
    # Filter out all genes with not enough guides pointing to them
    maude_gene_stats <- maude_gene_stats %>%
      filter(numGuides >= min_guides_per_gene)
    

    logger::log_info("Saving MAUDE gene stats to: {maude_gene_stats_fpath}")
    
    saveRDS(maude_gene_stats,maude_gene_stats_fpath)
  } else  if (isTRUE(run_plots_stage)){
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
    maude_bins = maude_bins,
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