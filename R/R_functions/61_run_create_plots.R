# 61_run_create_plots.R

run_create_plots <- function(opt, file_info_suffix) {
  
  logger::log_info("Starting plot generation...")
  
  logger::log_info("Loading Results...")
  loaded_results <- load_results_for_plotting(
    file_info_suffix = file_info_suffix,
    opt = opt
  )

  run_count_violin_plots(
    count_df_long = loaded_results$count_df_long,
    opt = opt,
    y_limit = opt$plots_violin_y_limit
  )

  run_maude_qc_plots(
    count_df_long = loaded_results$count_df_long,
    maude_counts_df = loaded_results$maude_counts_df,
    opt = opt,
    input_recovery = isTRUE(opt$recover_input)
  )
  
  run_waterfall_plot(
    opt = opt
  )
  logger::log_info("Finished plot generation.")
  
  invisible(TRUE)
}