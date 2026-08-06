# R/R_functions/EA_21_create_maude_qc_plot_objects.R

create_maude_qc_plot_objects <- function(count_df_long,
                                         maude_counts_df,
                                         cfg,
                                         input_recovery = FALSE) {
  
  all_bins <- cfg$bins$bin_name
  
  if (identical(cfg$counting$data_type, "umis")) {
    data_name <- "UMIs"
  } else if (identical(cfg$counting$data_type, "reads")) {
    data_name <- "reads"
  } else {
    stop_log("`counting.data_type` must be either 'umis' or 'reads'.")
  }
  
  if (isTRUE(input_recovery)) {
    title_suffix_1 <- "(without input recovery)"
    title_suffix_2 <- "(with input recovery)"
  } else {
    title_suffix_1 <- "(before MAUDE preparation)"
    title_suffix_2 <- "(ready for MAUDE)"
  }
  
  sample_sum_df <- count_df_long %>%
    dplyr::group_by(bin_name, exp) %>%
    dplyr::summarise(total_count = sum(count), .groups = "drop") %>%
    dplyr::mutate(
      Sample = interaction(bin_name, exp, drop = TRUE),
      Bin = bin_name
    )
  
  p_sum_per_sample <- ggplot2::ggplot(
    sample_sum_df,
    ggplot2::aes(x = Sample, y = total_count, fill = Bin)
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = paste("Sum of", data_name, "per sample"),
      x = "Sample",
      y = paste("Sum of", data_name)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
  
  sample_n_df <- count_df_long %>%
    dplyr::group_by(bin_name, exp) %>%
    dplyr::summarise(num_sgRNAs = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      Sample = interaction(bin_name, exp, drop = TRUE),
      Bin = bin_name
    )
  
  p_n_per_sample <- ggplot2::ggplot(
    sample_n_df,
    ggplot2::aes(x = Sample, y = num_sgRNAs, fill = Bin)
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = paste("Number of sgRNA entries per sample", title_suffix_1),
      x = "Sample",
      y = "Number of sgRNA entries"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
  
  bin_n_before_df <- count_df_long %>%
    dplyr::group_by(bin_name) %>%
    dplyr::summarise(num_sgRNAs = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(Bin = bin_name)
  
  p_n_per_bin_before <- ggplot2::ggplot(
    bin_n_before_df,
    ggplot2::aes(x = Bin, y = num_sgRNAs, fill = Bin)
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = paste("Number of sgRNA entries per bin", title_suffix_1),
      x = "Bin",
      y = "sgRNA entries"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
  
  bin_n_after_df <- maude_counts_df %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(all_bins),
        ~ sum(!is.na(.x))
      )
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "Bin",
      values_to = "count_non_na"
    )
  
  p_n_per_bin_after <- ggplot2::ggplot(
    bin_n_after_df,
    ggplot2::aes(x = Bin, y = count_non_na, fill = Bin)
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = paste("Number of sgRNA entries per bin", title_suffix_2),
      x = "Bin",
      y = "Number of sgRNA entries"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
  
  bin_sum_after_df <- maude_counts_df %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(all_bins),
        ~ sum(.x, na.rm = TRUE)
      )
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "Bin",
      values_to = "total_count"
    )
  
  p_sum_per_bin_after <- ggplot2::ggplot(
    bin_sum_after_df,
    ggplot2::aes(x = Bin, y = total_count, fill = Bin)
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = paste("Sum of", data_name, "per bin"),
      x = "Bin",
      y = "Total counts"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  manual_plots <- list(
    sum_per_sample = p_sum_per_sample,
    n_per_sample_before_maude = p_n_per_sample,
    n_per_bin_before_maude = p_n_per_bin_before,
    n_per_bin_ready_for_maude = p_n_per_bin_after,
    sum_per_bin_ready_for_maude = p_sum_per_bin_after
  )
  
  control_plots <- control_sanity_check(
    cfg = cfg,
    maude_counts_df = maude_counts_df,
    print = FALSE
  )
  
  maude_qc_plots <- plot_maude_qc(
    maude_counts_df = maude_counts_df,
    cfg = cfg, 
    print = FALSE
  )
  
  return(c(
    manual_plots,
    control_plots,
    maude_qc_plots
  ))
}