#-------------------------------------------------------------------------------
# Load function for count_df_long, maude_counts_df and hits_current_settings
#-------------------------------------------------------------------------------

load_results_for_plotting <- function(file_info_suffix, opt){
  #-----------------------------------------------------------------------------
  # Make count_df_long
  #-----------------------------------------------------------------------------
  logger::log_info("Reading in read/UMI counts...")
  
  if (identical(opt$read_counting, "bcwithqc")){
    count_df_long <- process_bcwithqc_data(
      data_type = opt$data_type,
      skip_list = opt$skip_list)
  } else if (identical(opt$read_counting, "align_UMI_tools")) {
    if (identical(opt$data_type, "umis")){
      count_df_long <- process_folder_files(dedup_output_folder,
                                            skip_list = skip_list) #Add threshold df if thresholds should be applied 
    }
    if (identical(opt$data_type, "reads")){
      count_df_long <- process_folder_files(mapped_output_folder,
                                            skip_list = opt$skip_list) #Add threshold df if thresholds should be applied 
    }
  } else {
    stop(opt$read_counting, " is not a valid read_counting value. Exiting. read_counting must be 'bcwithqc' or 'align_UMI_tools'.")
  }
  logger::log_info("Finished Reading in read/UMI counts.")
  logger::log_info("sgRNAs aligned to the wrong sublibrary were excluded from the analysis.")
  
  #-----------------------------------------------------------------------------
  # Optional count_df_long include_controls_list if use_only_these_controls_list is given
  #-----------------------------------------------------------------------------
  if (exists("use_only_these_controls_list")) {
    if (length(use_only_these_controls_list) > 0){
      # list of sgRNAs that were not targeting (before) and not in allowed controls
      excluded_controls <- count_df_long %>%
        filter(group_category != "targeting", !sgRNA %in% use_only_these_controls_list) %>%
        distinct(sgRNA) %>%
        pull(sgRNA)
      
      include_controls_list <- c(include_controls_list, excluded_controls)
    }
  }
  if (length(include_controls_list) > 0){
    for (control_gene in include_controls_list){
      count_df_long$group_category[grepl(control_gene, count_df_long$sgRNA)] <- "targeting"
    }
  }
  
  #-----------------------------------------------------------------------------
  # Optional Combine either samples or sublibraries of the same condition. 
  #-----------------------------------------------------------------------------
  
  if (combine_for_guide_stats != ""){
    if (combine_for_guide_stats == "sample"){
      count_df_long <- count_df_long %>%
        group_by(sgRNA, sublib, condition) %>%
        summarise(
          # set sample name for the combined rows
          sample = "sample_1",
          count = sum(count, na.rm = TRUE),
          exp = dplyr::first(exp),
          group_category = dplyr::first(group_category),
          
          # keep one sublib/sgRNA (also fine even though they're grouping keys)
          .groups = "drop"
        )
    }
    if (combine_for_guide_stats == "sublib"){
      count_df_long <- count_df_long %>%
        group_by(sgRNA, sample, condition) %>%
        summarise(
          sublib = "sublib_1",
          count = sum(count, na.rm = TRUE),
          exp = dplyr::first(exp),
          group_category = dplyr::first(group_category),
          
          # keep one sublib/sgRNA (also fine even though they're grouping keys)
          .groups = "drop"
        )    
    }
  }
  
  Hits_current_settings <- add_info_wrapper(file_info_suffix)
  if (length(include_controls_list) > 0) {
    for (control_gene in include_controls_list) {
      # here we remove everything after the last _, so stuff like AAVS1_9 and
      # AAVS1_13 are both treated as AAVS1
      control_gene <- sub("_[^_]*$", "", control_gene)
      Hits_current_settings$symbol[grepl(control_gene, Hits_current_settings$entrez)] <- control_gene
    }
    # Filter out rows where the 'symbol' is in the include_controls_list
    Hits_current_settings_no_controls <- Hits_current_settings %>% 
      filter(!symbol %in% include_controls_list | is.na(symbol))
  }
  
  maude_counts_df <- run_prepare_data_for_MAUDE(
    count_df_long = count_df_long,
    opt = opt)
  
  return(invisible(list(
    count_df_long = count_df_long,
    maude_counts_df = maude_counts_df,
    hits_current_settings = Hits_current_settings,
    hits_current_settings_no_controls = Hits_current_settings_no_controls
  )))
}

#-------------------------------------------------------------------------------
# post read counting Violin Plot generation
#-------------------------------------------------------------------------------

run_count_violin_plots <- function(count_df_long,
                                   opt,
                                   y_limit = 8000,
                                   non_targeting = TRUE,
                                   output_subdir = "01_read_or_umi_count_plots") {
  
  plot_dir <- file.path(plots_output_folder, output_subdir)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  logger::log_info("Creating read/UMI count violin plots in: {plot_dir}")
  
  # ---------------------------------------------------------------------------
  # Summary tables
  # ---------------------------------------------------------------------------
  
  summary_medians <- get_grouped_summary_wide(count_df_long, stat = "median")
  summary_means <- get_grouped_summary_wide(count_df_long, stat = "mean")
  
  summary_xlsx <- file.path(
    plot_dir,
    "count_summary.xlsx"
  )
  
  writexl::write_xlsx(
    list(
      medians = summary_medians,
      means = summary_means
    ),
    summary_xlsx
  )
  
  logger::log_info("Saved count summary tables to: {summary_xlsx}")
  
  # ---------------------------------------------------------------------------
  # Helper for saving named plot lists
  # ---------------------------------------------------------------------------
  
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
  
  # ---------------------------------------------------------------------------
  # Per sublib/sample violin plots
  # ---------------------------------------------------------------------------
  
  plots_raw <- plot_violin_by_sublib_sample(
    count_df_long,
    norm_method = NULL
  )
  
  save_plot_list(
    plot_list = plots_raw,
    prefix = "violin_by_sublib_sample_raw"
  )
  
  plots_norm <- plot_violin_by_sublib_sample(
    count_df_long,
    norm_method = "control_median"
  )
  
  save_plot_list(
    plot_list = plots_norm,
    prefix = "violin_by_sublib_sample_control_median"
  )
  
  # ---------------------------------------------------------------------------
  # Targeting / non-targeting split plots
  # ---------------------------------------------------------------------------
  
  targeting_raw <- plot_violin_by_group_category_split_by_sublib(
    count_df_long,
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
    count_df_long,
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
      count_df_long,
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
      count_df_long,
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
                               opt,
                               input_recovery = FALSE,
                               output_subdir = "02_MAUDE_QC_plots") {
  
  plot_dir <- file.path(plots_output_folder, output_subdir)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  logger::log_info("Creating MAUDE QC plots in: {plot_dir}")
  
  plots <- create_maude_qc_plot_objects(
    count_df_long = count_df_long,
    maude_counts_df = maude_counts_df,
    opt = opt,
    input_recovery = input_recovery
  )
  
  for (plot_name in names(plots)) {
    plot_obj <- plots[[plot_name]]
    
    if (!inherits(plot_obj, "ggplot")) {
      logger::log_warn("Skipping non-ggplot object in MAUDE QC plot list: {plot_name}")
      next
    }
    
    safe_name <- gsub("[^A-Za-z0-9_.-]+", "_", plot_name)
    
    out_path <- file.path(
      plot_dir,
      paste0(safe_name, ".png")
    )
    
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
create_waterfall_plot <- function(Hits_df, opt) {
  
  
  if (length(include_controls_list) > 0) {
    for (control_gene in include_controls_list) {
      # here we remove everything after the last _, so stuff like AAVS1_9 and
      # AAVS1_13 are both treated as AAVS1
      control_gene <- sub("_[^_]*$", "", control_gene)
      Hits_df$symbol[grepl(control_gene, Hits_df$entrez)] <- control_gene
    }
  }
  
  plot_df <- Hits_df
  
  p <- plot_significance_by_rank(
    plot_df,
    mark_cntrl = opt$plots_waterfall_mark_cntrl,
    mark_special = opt$plots_waterfall_mark_special,
    mark_N_top_hits = opt$plots_waterfall_mark_N_top_hits,
    box_padding = opt$plots_waterfall_box_padding,
    no_text = opt$plots_waterfall_no_text,
    signif_lines = opt$plots_waterfall_signif_lines,
    mark_all_signif_level = opt$plots_waterfall_mark_all_signif_level,
    break_in_plot = opt$plots_waterfall_break_in_plot,
    top_padding = opt$plots_waterfall_top_padding,
    custom_title = opt$plots_waterfall_custom_title
  )
  return(invisible(p))
}

run_waterfall_plot <- function(opt, output_subdir = "03_waterfall_plots") {
  
  plot_dir <- file.path(plots_output_folder, output_subdir)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  Hits_df <- add_info_wrapper(file_info_suffix)

  p <- create_waterfall_plot(Hits_df = Hits_df, opt = opt)
  
  out_path <- file.path(
    plot_dir,
    paste0(make_waterfall_filename(opt), ".", opt$plots_waterfall_file_format)
  )
  
  ggplot2::ggsave(
    filename = out_path,
    plot = p,
    width = opt$plots_waterfall_width,
    height = opt$plots_waterfall_height
  )
  
  logger::log_info("Saved waterfall plot to: {out_path}")
  
  invisible(p)
}
make_waterfall_filename <- function(opt) {
  
  parts <- c("waterfall")
  
  if (isTRUE(opt$plots_waterfall_mark_cntrl)) {
    parts <- c(parts, "controls")
  }
  
  mark_special <- opt$plots_waterfall_mark_special
  
  has_mark_special <- !is.null(mark_special) &&
    length(mark_special) > 0 &&
    !all(is.na(mark_special)) &&
    !all(mark_special == "")
  
  if (has_mark_special) {
    special <- paste(mark_special, collapse = "-")
    special <- gsub("[^A-Za-z0-9_.-]+", "_", special)
    parts <- c(parts, paste0("mark_", special))
  }
  
  if (opt$plots_waterfall_mark_N_top_hits > 0) {
    parts <- c(parts, paste0("top", opt$plots_waterfall_mark_N_top_hits))
  }
  
  if (!is.null(opt$plots_waterfall_mark_all_signif_level) &&
      !is.na(opt$plots_waterfall_mark_all_signif_level)) {
    parts <- c(parts, paste0("FDR", opt$plots_waterfall_mark_all_signif_level))
  }
  
  if (isTRUE(opt$plots_waterfall_no_text)) {
    parts <- c(parts, "no_text")
  }
  
  return(invisible(paste(parts, collapse = "_")))
}