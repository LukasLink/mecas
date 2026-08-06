# R/R_functions/EA_30_run_waterfall_plot.R

#-------------------------------------------------------------------------------
# post MAUDE Waterfall plot
#-------------------------------------------------------------------------------

run_waterfall_plot <- function(cfg,
                               output_subdir = "03_waterfall_plots") {
  
  plot_dir <- file.path(cfg$paths$plots_output_folder, output_subdir)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  Hits_df <- add_info_wrapper(
    file_suffix = cfg$suffix$file_suffix,
    cfg = cfg
  )
  
  p <- create_waterfall_plot(
    Hits_df = Hits_df,
    cfg = cfg
  )
  
  out_path <- file.path(
    plot_dir,
    paste0(
      make_waterfall_filename(cfg),
      ".",
      cfg$plots$waterfall$file_format
    )
  )
  
  ggplot2::ggsave(
    filename = out_path,
    plot = p,
    width = cfg$plots$waterfall$width,
    height = cfg$plots$waterfall$height
  )
  
  logger::log_info("Saved waterfall plot to: {out_path}")
  
  invisible(p)
}
