# run_plot.R

run_plot <- function(cfg) {
  
  logger::log_info("Starting plot generation...")
  
  # Quick Hack for nicer output names
  cfg$suffix$file_suffix <- ".rds"
  cfg$suffix$file_info_suffix <- ""
  
  logger::log_info("Loading Results...")
  
  loaded_results <- load_results_for_plotting(
    file_suffix = cfg$suffix$file_suffix,
    cfg = cfg
  )
  
  run_count_violin_plots(
    count_df_long = loaded_results$count_df_long,
    cfg = cfg,
    y_limit = cfg$plots$read_count_violin$violin_y_limit
  )
  
  run_maude_qc_plots(
    count_df_long = loaded_results$count_df_long,
    maude_counts_df = loaded_results$maude_counts_df,
    cfg = cfg,
    input_recovery = isTRUE(cfg$normalization$recover_input)
  )
  
  run_waterfall_plot(
    cfg = cfg
  )
  
  logger::log_info("Finished plot generation.")
  
  invisible(TRUE)
}