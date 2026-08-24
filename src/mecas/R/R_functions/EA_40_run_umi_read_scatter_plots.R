# R/R_functions/EA_40_run_umi_read_scatter_plots.R

run_umi_read_scatter_plots <- function(umi_counts_df, read_counts_df, cfg){
  
  plot_dir <- file.path(cfg$paths$plots_output_folder, "01_read_or_umi_count_plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  umi_counts_df <- umi_counts_df %>% rename(umi_count = count)
  read_counts_df <- read_counts_df %>% rename(read_count = count)

  # Join read + UMI counts and add gene
  combined_df <- full_join(
    read_counts_df,
    umi_counts_df,
    by = c("sgRNA", "sample", "bin_name", "sublib")
    ) %>%
    mutate(
      read_count = replace_na(read_count, 0),
      umi_count  = replace_na(umi_count, 0)
    ) %>% 
    mutate(
      read_count = read_count + 1,
      umi_count = umi_count + 1)
  
  plot_specs <- combined_df %>%
    group_by(sample, bin_name, sublib) %>%
    group_split()
  
  for (df in plot_specs) {
    spec <- scatter_plot_for_one_group(df)
    
    ggsave(
      filename = file.path(plot_dir, spec$filename),
      plot = spec$plot,
      width = 6,
      height = 5,
      dpi = 300
    )
  }
  
  logger::log_info("Saved UMI/read scatter plots to: {plot_dir}")
  
  invisible(list(
    plot_dir = plot_dir,
    n_plots = length(plot_specs)
  ))
}