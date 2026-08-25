# R/R_functions/EA_60_run_reads_per_umi_distribution_plots.R

run_reads_per_umi_distribution_plots <- function(count_df_long,
                                                 cfg,
                                                 display_quantile = 0.99,
                                                 output_subdir = "01a_reads_per_umi_distribution_plots") {
  
  plot_dir <- file.path(cfg$paths$plots_output_folder, output_subdir)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  logger::log_info("Creating reads-per-UMI distribution plots in: {plot_dir}")
  
  # Collapse the potentially very large dataframe before plotting.
  # Each resulting row represents one observed reads-per-UMI count within a sample/bin.
  data.table::setDT(count_df_long)
  
  distribution_df <- count_df_long[
    ,
    .(n_umis = .N),
    by = .(sample, bin_name, count)
  ]
  
  data.table::setorder(distribution_df, sample, bin_name, count)
  
  distribution_df[
    ,
    `:=`(
      fraction = n_umis / sum(n_umis),
      cumulative_fraction = cumsum(n_umis) / sum(n_umis)
    ),
    by = .(sample, bin_name)
  ]
  
  # Determine the requested upper display limit separately for each sample/bin.
  quantile_limits <- distribution_df[
    cumulative_fraction >= display_quantile,
    .SD[1],
    by = .(sample, bin_name)
  ][
    ,
    .(sample, bin_name, display_limit = count)
  ]
  
  distribution_df <- merge(
    distribution_df,
    quantile_limits,
    by = c("sample", "bin_name"),
    all.x = TRUE
  )
  
  distribution_99_df <- distribution_df[count <= display_limit]
  
  # ---------------------------------------------------------------------------
  # Reads-per-UMI frequency distribution
  # ---------------------------------------------------------------------------
  
  p_distribution <- ggplot2::ggplot(
    distribution_99_df,
    ggplot2::aes(x = factor(count), y = fraction, fill = bin_name)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~ sample, scales = "free_x") +
    ggplot2::labs(
      title = paste0("Reads per UMI distribution (first ", display_quantile * 100, "% of UMIs)"),
      x = "Reads per UMI",
      y = "Fraction of bins",
      color = "Bin"
    ) +
    ggplot2::theme_bw()
  
  distribution_fpath <- file.path(plot_dir, "reads_per_umi_distribution_99_percent.png")
  
  ggplot2::ggsave(
    filename = distribution_fpath,
    plot = p_distribution,
    width = 9,
    height = 6,
    dpi = 300
  )
  
  # ---------------------------------------------------------------------------
  # Reads-per-UMI ECDF
  # ---------------------------------------------------------------------------
  
  p_ecdf <- ggplot2::ggplot(
    distribution_df,
    ggplot2::aes(x = count, y = cumulative_fraction, color = bin_name, group = bin_name)
  ) +
    ggplot2::geom_step(linewidth = 0.7) +
    ggplot2::facet_wrap(~ sample, scales = "free_x") +
    ggplot2::labs(
      title = "Cumulative reads per UMI distribution",
      x = "Reads per UMI",
      y = "Cumulative fraction of UMIs",
      color = "Bin"
    ) +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(limits = c(0, 1))
  
  ecdf_fpath <- file.path(plot_dir, "reads_per_umi_ecdf.png")
  
  ggplot2::ggsave(
    filename = ecdf_fpath,
    plot = p_ecdf,
    width = 9,
    height = 6,
    dpi = 300
  )
  
  # ---------------------------------------------------------------------------
  # Export underlying distribution table
  # ---------------------------------------------------------------------------
  
  distribution_table_fpath <- file.path(plot_dir, "reads_per_umi_distribution.tsv.gz")
  
  readr::write_tsv(
    as.data.frame(distribution_df),
    distribution_table_fpath
  )
  
  logger::log_info("Saved reads-per-UMI distribution plot to: {distribution_fpath}")
  logger::log_info("Saved reads-per-UMI ECDF plot to: {ecdf_fpath}")
  logger::log_info("Saved reads-per-UMI distribution table to: {distribution_table_fpath}")
  logger::log_info("Finished creating reads-per-UMI distribution plots.")
  
  invisible(list(
    plot_dir = plot_dir,
    distribution_plot = distribution_fpath,
    ecdf_plot = ecdf_fpath,
    distribution_table = distribution_table_fpath
  ))
}