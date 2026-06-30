run_consensus_call <- function(maude_counts_df,
                               cfg,
                               run_maude_stage,
                               run_plots_stage) {

  
  n_reps <- cfg$consensus$n_reps
  high_confidence_fdr <- cfg$consensus$high_confidence$FDR_threshold
  high_confidence_hits_in_reps <- cfg$consensus$high_confidence$hits_in_X_reps
  explorative_fdr <- cfg$consensus$explorative$FDR_threshold
  explorative_hits_in_reps <- cfg$consensus$explorative$hits_in_X_reps
  
  logger::log_info("Starting consensus MAUDE calling with {n_reps} replicate runs.")
  
  base_file_suffix <- cfg$suffix$file_suffix
  
  for (rep_i in seq_len(n_reps)) {
    
    rep_suffix <- sub(
      "\\.rds$",
      paste0("_rep", rep_i, ".rds"),
      base_file_suffix
    )
    
    logger::log_info("Running MAUDE replicate {rep_i}/{n_reps} with suffix: {rep_suffix}")
    
    run_MAUDE(
      maude_counts_df = maude_counts_df,
      cfg = cfg,
      file_suffix = rep_suffix,
      run_maude_stage = run_maude_stage,
      run_plots_stage = run_plots_stage
    )
  }
  
  logger::log_info("Finished MAUDE replicate runs. Starting high-confidence consensus call.")
  
  high_confidence_hits <- automate_calc_replicate_comparison(
    cfg = cfg,
    FDR_threshold = high_confidence_fdr,
    hits_in_X_reps = high_confidence_hits_in_reps,
    correlation_heatmap = FALSE,
    venn_diagram = FALSE
  )
  
  logger::log_info("Starting explorative consensus call.")
  
  explorative_hits <- automate_calc_replicate_comparison(
    cfg = cfg,
    FDR_threshold = explorative_fdr,
    hits_in_X_reps = explorative_hits_in_reps,
    correlation_heatmap = FALSE,
    venn_diagram = FALSE
  )
  
  if (cfg$suffix$file_info_suffix == ""){
    high_confidence_xlsx <- get_file_path(
      cfg$paths$results_output_folder,
      paste0("Hits_high_confidence", cfg$suffix$file_info_suffix, ".xlsx")
    )
    
    high_confidence_rds <- get_file_path(
      cfg$paths$rds_output_folder,
      paste0("high_confidence_hits", cfg$suffix$file_info_suffix, ".rds")
    )
    
    explorative_xlsx <- get_file_path(
      cfg$paths$results_output_folder,
      paste0("Hits_explorative", cfg$suffix$file_info_suffix, ".xlsx")
    )
    
    explorative_rds <- get_file_path(
      cfg$paths$rds_output_folder,
      paste0("explorative_hits", cfg$suffix$file_info_suffix, ".rds")
    )
  } else {
    high_confidence_xlsx <- get_file_path(
      cfg$paths$results_output_folder,
      paste0("Hits_high_confidence_", cfg$suffix$file_info_suffix, ".xlsx")
    )
    
    high_confidence_rds <- get_file_path(
      cfg$paths$rds_output_folder,
      paste0("high_confidence_hits_", cfg$suffix$file_info_suffix, ".rds")
    )
    
    explorative_xlsx <- get_file_path(
      cfg$paths$results_output_folder,
      paste0("Hits_explorative_", cfg$suffix$file_info_suffix, ".xlsx")
    )
    
    explorative_rds <- get_file_path(
      cfg$paths$rds_output_folder,
      paste0("explorative_hits_", cfg$suffix$file_info_suffix, ".rds")
    )
  }
  

  
  writexl::write_xlsx(high_confidence_hits$Hits_in_X_df, high_confidence_xlsx)
  saveRDS(high_confidence_hits$Hits_in_X_df, high_confidence_rds)
  
  writexl::write_xlsx(explorative_hits$Hits_in_X_df, explorative_xlsx)
  saveRDS(explorative_hits$Hits_in_X_df, explorative_rds)
  
  logger::log_info("Saved high-confidence hits to: {high_confidence_xlsx}")
  logger::log_info("Saved explorative hits to: {explorative_xlsx}")
  
  logger::log_info("Calculating replicate means.")
  
  replicate_means <- automate_calc_replicate_means(cfg = cfg)
  
  logger::log_info("Finished consensus calling.")
  
  invisible(list(
    high_confidence_hits = high_confidence_hits,
    explorative_hits = explorative_hits,
    replicate_means = replicate_means,
    high_confidence_xlsx = high_confidence_xlsx,
    explorative_xlsx = explorative_xlsx,
    high_confidence_rds = high_confidence_rds,
    explorative_rds = explorative_rds
  ))
}