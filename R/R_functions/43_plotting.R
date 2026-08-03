#-------------------------------------------------------------------------------
# post read counting Violin Plot generation
#-------------------------------------------------------------------------------

run_count_violin_plots <- function(count_df_long,
                                   cfg,
                                   y_limit = 8000,
                                   non_targeting = TRUE,
                                   output_subdir = "01_read_or_umi_count_plots") {
  
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
    prefix = "violin_by_sublib_sample_raw"
  )
  
  plots_norm <- plot_violin_by_sublib_sample(
    count_df_long = count_df_long,
    cfg = cfg,
    norm_method = "control_median"
  )
  
  save_plot_list(
    plot_list = plots_norm,
    prefix = "violin_by_sublib_sample_control_median"
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

#-------------------------------------------------------------------------------
# post MAUDE Waterfall plot
#-------------------------------------------------------------------------------

create_waterfall_plot <- function(Hits_df, cfg) {
  
  include_controls_list <- cfg$controls$include_controls %||% character()
  
  if (length(include_controls_list) > 0) {
    for (control_gene in include_controls_list) {
      control_gene <- sub("_[^_]*$", "", control_gene)
      Hits_df$symbol[
        grepl(control_gene, Hits_df$entrez)
      ] <- control_gene
    }
  }
  
  plot_df <- Hits_df
  
  p <- plot_significance_by_rank(
    plot_df,
    mark_cntrl = cfg$plots$waterfall$mark_cntrl,
    mark_special = cfg$plots$waterfall$mark_special,
    mark_N_top_hits = cfg$plots$waterfall$mark_N_top_hits,
    box_padding = cfg$plots$waterfall$box_padding,
    no_text = cfg$plots$waterfall$no_text,
    signif_lines = cfg$plots$waterfall$signif_lines,
    mark_all_signif_level = cfg$plots$waterfall$mark_all_signif_level,
    break_in_plot = cfg$plots$waterfall$break_in_plot,
    top_padding = cfg$plots$waterfall$top_padding,
    custom_title = cfg$plots$waterfall$custom_title
  )
  
  invisible(p)
}

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

make_waterfall_filename <- function(cfg) {
  
  parts <- c("waterfall")
  
  if (isTRUE(cfg$plots$waterfall$mark_cntrl)) {
    parts <- c(parts, "controls")
  }
  
  mark_special <- cfg$plots$waterfall$mark_special
  
  has_mark_special <- !is.null(mark_special) &&
    length(mark_special) > 0 &&
    !all(is.na(mark_special)) &&
    !all(mark_special == "")
  
  if (has_mark_special) {
    special <- paste(mark_special, collapse = "-")
    special <- gsub("[^A-Za-z0-9_.-]+", "_", special)
    parts <- c(parts, paste0("mark_", special))
  }
  
  if (cfg$plots$waterfall$mark_N_top_hits > 0) {
    parts <- c(parts, paste0("top", cfg$plots$waterfall$mark_N_top_hits))
  }
  
  if (!is.null(cfg$plots$waterfall$mark_all_signif_level) &&
      !is.na(cfg$plots$waterfall$mark_all_signif_level)) {
    parts <- c(parts, paste0("FDR", cfg$plots$waterfall$mark_all_signif_level))
  }
  
  if (isTRUE(cfg$plots$waterfall$no_text)) {
    parts <- c(parts, "no_text")
  }
  
  invisible(paste(parts, collapse = "_"))
}