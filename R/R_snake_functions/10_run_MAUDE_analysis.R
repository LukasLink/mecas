# R/R_functions_snake/10_run_MAUDE.R

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

  # Prepare data for MAUDE
  maude_counts_df <- run_prepare_data_for_MAUDE(
    count_df_long = count_df_long,
    cfg = cfg
  )
  
  logger::log_info("Starting initial MAUDE run ...")
  
  maude_results <- run_MAUDE(
    maude_counts_df = maude_counts_df,
    cfg = cfg,
    file_suffix = cfg$suffix$file_suffix,
    run_maude_stage = TRUE,
    run_plots_stage = TRUE
  )
  
  # Export MAUDE results + combine replicates
  export_df <- handle_auto_combine_replicates_and_export(
    file_suffix = cfg$suffix$file_suffix,
    maude_results = maude_results,
    cfg = cfg
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

  logger::log_info("DONE!")
  
  invisible(cfg)
}