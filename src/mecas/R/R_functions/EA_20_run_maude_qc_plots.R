# R/R_functions/EA_20_run_maude_qc_plots.R

#-------------------------------------------------------------------------------
# pre MAUDE QC plots generation
#-------------------------------------------------------------------------------

run_maude_qc_plots <- function(count_df_long,
                               maude_counts_df,
                               cfg,
                               input_recovery = FALSE,
                               output_subdir = "02_MAUDE_QC_plots") {
  
  plot_dir <- file.path(cfg$paths$plots_output_folder, output_subdir)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  logger::log_info("Creating MAUDE QC plots in: {plot_dir}")
  
  plots <- create_maude_qc_plot_objects(
    count_df_long = count_df_long,
    maude_counts_df = maude_counts_df,
    cfg = cfg,
    input_recovery = input_recovery
  )
  
  for (plot_name in names(plots)) {
    plot_obj <- plots[[plot_name]]
    
    if (!inherits(plot_obj, "ggplot")) {
      logger::log_warn("Skipping non-ggplot object in MAUDE QC plot list: {plot_name}")
      next
    }
    
    safe_name <- gsub("[^A-Za-z0-9_.-]+", "_", plot_name)
    
    out_path <- file.path(plot_dir, paste0(safe_name, ".png"))
    
    ggplot2::ggsave(
      filename = out_path,
      plot = plot_obj,
      width = 8,
      height = 5,
      dpi = 300
    )
  }
  
  logger::log_info("Finished MAUDE QC plots.")
  
  invisible(plots)
}