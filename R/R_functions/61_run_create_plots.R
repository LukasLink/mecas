# 61_run_create_plots.R
# THIS IS DEPRACTED AND REPLACED WITH run_plot
run_create_plots <- function(cfg) {
  
  logger::log_info("Starting plot generation...")
  
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