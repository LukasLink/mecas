# run_plot.R

run_plot <- function(cfg) {
  
  logger::log_info("Starting plot generation...")
  
  # Quick Hack for nicer output names
  cfg$suffix$file_suffix <- ".rds"
  cfg$suffix$file_info_suffix <- ""
  
  logger::log_info("Loading Results...")

  count_df_long <- readRDS(cfg$paths$count_df_fpath)
  
  maude_counts_df <- readRDS(cfg$paths$maude_counts_df_fpath)
  
  run_count_violin_plots(
    count_df_long = count_df_long,
    cfg = cfg,
    y_limit = cfg$plots$read_count_violin$violin_y_limit
  )
  
  run_maude_qc_plots(
    count_df_long = count_df_long,
    maude_counts_df = maude_counts_df,
    cfg = cfg,
    input_recovery = isTRUE(cfg$normalization$recover_input)
  )
  
  run_waterfall_plot(
    cfg = cfg
  )
  
  logger::log_info("Finished plot generation.")
  
  invisible(TRUE)
}