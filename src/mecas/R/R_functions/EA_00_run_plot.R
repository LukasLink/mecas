# R/R_functions/EA_00_run_plot.R

run_plot <- function(cfg) {
  
  logger::log_info("Starting plot generation...")
  
  # Quick Hack for nicer output names
  cfg$suffix$file_suffix <- ".rds"
  cfg$suffix$file_info_suffix <- ""
  
  if (isTRUE(cfg$counting$umis_as_sublibs)) {
    run_count_violin_plots_umis_as_sublibs(
      count_df_long = readRDS(cfg$paths$count_df_fpath),
      cfg = cfg,
      y_limit = cfg$plots$read_count_violin$violin_y_limit
    )
    run_reads_per_umi_distribution_plots(
      count_df_long = readRDS(cfg$paths$count_df_fpath),
      cfg = cfg
    )
  } else {
    run_count_violin_plots(
      count_df_long = readRDS(cfg$paths$count_df_fpath),
      cfg = cfg,
      y_limit = cfg$plots$read_count_violin$violin_y_limit
    )
  }
  
  if (file.exists(cfg$paths$maude_counts_df_fpath)){
    
    run_maude_qc_plots(
      count_df_long = readRDS(cfg$paths$count_df_fpath),
      maude_counts_df = readRDS(cfg$paths$maude_counts_df_fpath),
      cfg = cfg,
      input_recovery = isTRUE(cfg$normalization$recover_input)
    )
    
    run_waterfall_plot(
      cfg = cfg
    )
  } else {
    log_info("No MAUDE results found at:\n{cfg$paths$maude_counts_df_fpath}\n Skipping MAUDE QC and Waterfall Hits plot generation.")
  }

  log_info("We are making duplicated plots in run_umi_read_scatter")
  if (identical("umis", cfg$counting$data_type) & file.exists(cfg$paths$reads_count_df_fpath)){
    if (isTRUE(cfg$counting$umis_as_sublibs)){
      run_umi_read_scatter_plots(
        umi_counts_df = get_umi_counts_in_umis_as_sublibs_mode(cfg = cfg),
        read_counts_df = readRDS(cfg$paths$reads_count_df_fpath),
        cfg = cfg
      )
    } else {
      run_umi_read_scatter_plots(
        umi_counts_df = readRDS(cfg$paths$count_df_fpath),
        read_counts_df = readRDS(cfg$paths$reads_count_df_fpath),
        cfg = cfg
      )
    }
  }
  log_info("We are making duplicated plots in run_count_density")
  if (identical("umis", cfg$counting$data_type)){
    run_count_density_plots(
      count_df_long = get_umi_counts_in_umis_as_sublibs_mode(cfg = cfg),
      cfg = cfg
    )
  } else {
    run_count_density_plots(
      count_df_long = readRDS(cfg$paths$count_df_fpath),
      cfg = cfg
    )
  }
  
  safe_copy_bcwithqc_plots(cfg = cfg)
  
  logger::log_info("Finished plot generation.")
  
  invisible(TRUE)
}