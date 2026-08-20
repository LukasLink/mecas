# R/R_functions/DA_00_run_MAUDE_analyis.R

run_MAUDE_analysis <- function(cfg){
  
  count_df_long <- tryCatch(
    {readRDS(cfg$paths$count_df_fpath)},
    error = function(e) {
      stop_log(
        paste0(
          "Failed to load count dataframe from: ",
          cfg$paths$count_df_fpath,
          "\nOriginal error: ",
          conditionMessage(e)
        )
      )
    }
  )
  
  # Quick Hack for nicer output names
  cfg$suffix$file_suffix <- ".rds"
  cfg$suffix$file_info_suffix <- ""
  
  #-----------------------------------------------------------------------------
  # Run MAUDE
  #-----------------------------------------------------------------------------
  log_info(
    "Before run_prepare_data_for_MAUDE | R memory: {sprintf('%.2f GB', memory_used_gb())} | RSS: {sprintf('%.2f GB', memory_rss_gb())}"
  )
  # Prepare data for MAUDE
  maude_counts_df <- run_prepare_data_for_MAUDE(
    count_df_long = count_df_long,
    cfg = cfg
  )
  log_info(
    "After run_prepare_data_for_MAUDE | R memory: {sprintf('%.2f GB', memory_used_gb())} | RSS: {sprintf('%.2f GB', memory_rss_gb())}"
  )
  logger::log_info("Starting initial MAUDE run ...")
  
  maude_results <- run_MAUDE(
    maude_counts_df = maude_counts_df,
    cfg = cfg,
    file_suffix = cfg$suffix$file_suffix,
    run_maude_stage = TRUE,
    run_plots_stage = FALSE
  )
  log_info(
    "After run_MAUDE | R memory: {sprintf('%.2f GB', memory_used_gb())} | RSS: {sprintf('%.2f GB', memory_rss_gb())}"
  )
  # Export MAUDE results + combine replicates
  export_df <- handle_auto_combine_replicates_and_export(
    file_suffix = cfg$suffix$file_suffix,
    maude_results = maude_results,
    cfg = cfg
  )
  log_info(
    "After export | R memory: {sprintf('%.2f GB', memory_used_gb())} | RSS: {sprintf('%.2f GB', memory_rss_gb())}"
  )  
  #-----------------------------------------------------------------------------
  # Finding Consensus Hits: MAUDE
  #-----------------------------------------------------------------------------
    
  if (isTRUE(cfg$consensus$run)) {
    consensus_results <- run_consensus_call(
      maude_counts_df = maude_counts_df,
      cfg = cfg,
      run_maude_stage = TRUE,
      run_plots_stage = TRUE
    )
  }
  log_info(
    "After run_consensus call | R memory: {sprintf('%.2f GB', memory_used_gb())} | RSS: {sprintf('%.2f GB', memory_rss_gb())}"
  )
  logger::log_info("DONE!")
  
  invisible(cfg)
}