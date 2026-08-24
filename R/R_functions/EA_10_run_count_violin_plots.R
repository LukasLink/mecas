# R/R_functions/EA_10_run_count_violin_plots.R

#-------------------------------------------------------------------------------
# post read counting Violin Plot generation
#-------------------------------------------------------------------------------

run_count_violin_plots <- function(count_df_long,
                                   cfg,
                                   y_limit = NA,
                                   auto_y_limit_prob =1,
                                   non_targeting = TRUE,
                                   output_subdir = "01_read_or_umi_count_plots") {
  
  get_auto_y_limit <- function(df, prob = auto_y_limit_prob) {
    if (!"count" %in% colnames(df)) {
      stop("count_df_long must contain a `count` column.")
    }
    
    y <- df$count
    y <- y[is.finite(y)]
    
    if (length(y) == 0) {
      stop("No finite count values found for automatic y_limit calculation.")
    }
    
    ceiling(stats::quantile(y, probs = prob, na.rm = TRUE, names = FALSE))
  }
  
  if (is.na(y_limit)) {
    y_limit <- get_auto_y_limit(count_df_long, prob = auto_y_limit_prob)
    logger::log_info("Auto-determined y_limit from data: {y_limit}")
  }
  
  plot_dir <- file.path(cfg$paths$plots_output_folder, output_subdir)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  logger::log_info("Creating read/UMI count violin plots in: {plot_dir}")
  
  summary_medians <- get_grouped_summary_wide(count_df_long, stat = "median")
  summary_means <- get_grouped_summary_wide(count_df_long, stat = "mean")
  
  summary_xlsx <- file.path(plot_dir, "count_summary.xlsx")
  
  writexl::write_xlsx(
    list(
      medians = summary_medians,
      means = summary_means
    ),
    summary_xlsx
  )
  
  logger::log_info("Saved count summary tables to: {summary_xlsx}")
  
  save_plot_list <- function(plot_list, prefix, width = 8, height = 5) {
    if (length(plot_list) == 0) {
      logger::log_warn("No plots generated for prefix: {prefix}")
      return(invisible(NULL))
    }
    
    for (plot_name in names(plot_list)) {
      safe_name <- gsub("[^A-Za-z0-9_.-]+", "_", plot_name)
      
      out_path <- file.path(
        plot_dir,
        paste0(prefix, "_", safe_name, ".png")
      )
      
      ggplot2::ggsave(
        filename = out_path,
        plot = plot_list[[plot_name]],
        width = width,
        height = height,
        dpi = 300
      )
    }
    
    invisible(NULL)
  }
  
  plots_raw <- plot_violin_by_sublib_sample(
    count_df_long = count_df_long,
    cfg = cfg,
    norm_method = NULL
  )
  
  save_plot_list(
    plot_list = plots_raw,
    prefix = "violin_counts_by_sublib_sample_raw"
  )
  
  plots_norm <- plot_violin_by_sublib_sample(
    count_df_long = count_df_long,
    cfg = cfg,
    norm_method = "control_median"
  )
  
  save_plot_list(
    plot_list = plots_norm,
    prefix = "violin_counts_by_sublib_sample_normalized"
  )
  
  targeting_raw <- plot_violin_by_group_category_split_by_sublib(
    df = count_df_long,
    cfg = cfg,
    include_targeting = TRUE,
    norm_method = NULL,
    y_limit = y_limit,
    viol_col = "lightgreen"
  )
  
  save_plot_list(
    plot_list = targeting_raw$plots,
    prefix = "targeting_raw",
    width = 10,
    height = 6
  )
  
  targeting_norm <- plot_violin_by_group_category_split_by_sublib(
    df = count_df_long,
    cfg = cfg,
    include_targeting = TRUE,
    norm_method = "control_median",
    y_limit = y_limit,
    viol_col = "lightgreen"
  )
  
  save_plot_list(
    plot_list = targeting_norm$plots,
    prefix = "targeting_control_median",
    width = 10,
    height = 6
  )
  
  if (isTRUE(non_targeting)) {
    
    non_targeting_raw <- plot_violin_by_group_category_split_by_sublib(
      df = count_df_long,
      cfg = cfg,
      include_targeting = FALSE,
      norm_method = NULL,
      y_limit = y_limit,
      viol_col = "lightblue"
    )
    
    save_plot_list(
      plot_list = non_targeting_raw$plots,
      prefix = "non_targeting_raw",
      width = 10,
      height = 6
    )
    
    non_targeting_norm <- plot_violin_by_group_category_split_by_sublib(
      df = count_df_long,
      cfg = cfg,
      include_targeting = FALSE,
      norm_method = "control_median",
      y_limit = y_limit,
      viol_col = "lightblue"
    )
    
    save_plot_list(
      plot_list = non_targeting_norm$plots,
      prefix = "non_targeting_control_median",
      width = 10,
      height = 6
    )
  }
  
  logger::log_info("Finished creating read/UMI count violin plots.")
  
  invisible(list(
    plot_dir = plot_dir,
    summary_xlsx = summary_xlsx
  ))
}
