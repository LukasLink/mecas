# R/R_functions/EA_40_run_count_density_plots.R

run_count_density_plots <- function(count_df_long, cfg, output_subdir = "01_read_or_umi_count_plots") {

  plot_dir <- file.path(cfg$paths$plots_output_folder, output_subdir)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  make_safe_name <- function(x) {
    x <- ifelse(is.na(x) | x == "", "NA", x)
    gsub("[^A-Za-z0-9_.-]+", "-", x)
  }
  
  # Build the plotting data
  if (identical("umis", cfg$counting$data_type)) {
    umi_counts_df <- count_df_long %>%
      mutate(count = replace_na(count, 0)) %>% 
      mutate(count = count + 1) %>% 
      rename(umi_count = count)
    
    read_counts_df <- readRDS(cfg$paths$reads_count_df_fpath) %>%
      mutate(count = replace_na(count, 0)) %>% 
      mutate(count = count + 1) %>% 
      rename(read_count = count)
    
    combined_df <- full_join(
      read_counts_df,
      umi_counts_df,
      by = c("sgRNA", "sample", "bin_name", "sublib")
    )
    
    plot_df <- combined_df %>%
      pivot_longer(
        cols = c(read_count, umi_count),
        names_to = "count_type",
        values_to = "count"
      ) %>%
      mutate(
        count_type = recode(
          count_type,
          read_count = "Read count",
          umi_count  = "UMI count"
        )
      )
  } else {
    plot_df <- count_df_long %>%
      mutate(count = replace_na(count, 0)) %>% 
      mutate(count = count + 1) %>% 
      mutate(count_type = "Count")
  }
  
  plot_groups <- plot_df %>%
    group_by(sample, bin_name, sublib) %>%
    group_split()
  
  for (df in plot_groups) {
    sample_name <- unique(df$sample)[1]
    bin_name <- unique(df$bin_name)[1]
    sublib <- unique(df$sublib)[1]
    
    title_txt <- paste(
      "Count distributions",
      paste0("sample = ", sample_name),
      paste0("bin_name = ", bin_name),
      paste0("sublib = ", sublib),
      sep = "\n"
    )
    
    if (identical("umis", cfg$counting$data_type)) {
      p <- ggplot(df, aes(x = count, fill = count_type)) +
        geom_density(alpha = 0.35) +
        scale_x_log10() +
        labs(
          title = title_txt,
          x = "Count + 1 (log10 scale)",
          y = "Density",
          fill = ""
        ) +
        theme_bw()
    } else {
      p <- ggplot(df, aes(x = count)) +
        geom_density(fill = "grey70", alpha = 0.6) +
        scale_x_log10() +
        labs(
          title = title_txt,
          x = "Count + 1 (log10 scale)",
          y = "Density"
        ) +
        theme_bw()
    }
    
    fname <- paste0(
      "count_density_",
      make_safe_name(sample_name), "_",
      make_safe_name(bin_name), "_",
      make_safe_name(sublib),
      ".png"
    )
    
    ggsave(
      filename = file.path(plot_dir, fname),
      plot = p,
      width = 6,
      height = 5,
      dpi = 300
    )
  }
  
  logger::log_info("Saved count density plots to: {plot_dir}")
  
  invisible(list(
    plot_dir = plot_dir,
    n_plots = length(plot_groups)
  ))
}