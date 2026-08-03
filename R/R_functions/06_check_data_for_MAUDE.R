plot_maude_qc <- function(maude_counts_df, cfg, print = TRUE, add_0s = FALSE) {
  
  all_bins <- cfg$bins$bin_name
  
  sorted_bins <- cfg$bins %>%
    dplyr::filter(sorted_or_unsorted == "sorted") %>%
    dplyr::pull(bin_name)
  
  if (isTRUE(add_0s)) {
    maude_counts_df <- maude_counts_df %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(sorted_bins),
          ~ ifelse(is.na(.x), 0, .x)
        )
      )
  }
  
  na_summary_overall <- maude_counts_df %>%
    pivot_longer(
      cols = dplyr::all_of(all_bins),
      names_to = "bin_name",
      values_to = "value"
    ) %>%
    mutate(is_na = is.na(value)) %>%
    group_by(bin_name) %>%
    summarise(na_count = sum(is_na), .groups = "drop")
  
  p1 <- ggplot(na_summary_overall, aes(x = bin_name, y = na_count, fill = bin_name)) +
    geom_col() +
    labs(title = "Missing Entries per Bin", y = "Number of NAs", x = "Bin") +
    theme_bw()
  
  na_summary_by_exp <- maude_counts_df %>%
    pivot_longer(
      cols = dplyr::all_of(all_bins),
      names_to = "bin_name",
      values_to = "value"
    ) %>%
    mutate(is_na = is.na(value)) %>%
    group_by(exp, bin_name) %>%
    summarise(na_count = sum(is_na), .groups = "drop")
  
  p2 <- ggplot(na_summary_by_exp, aes(x = exp, y = na_count, fill = bin_name)) +
    geom_col(position = "dodge") +
    labs(title = "Missing Entries per Bin by Experiment", x = "Experiment", y = "Number of NAs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  sgRNA_counts_per_exp <- maude_counts_df %>%
    group_by(exp, sgRNA) %>%
    summarise(dup_count = n(), .groups = "drop")
  
  p3 <- ggplot(sgRNA_counts_per_exp, aes(x = exp, y = dup_count)) +
    geom_boxplot(fill = "skyblue") +
    labs(
      title = "How Often Do Guides Appear Per Sample?",
      x = "Experiment",
      y = "Duplicated Guides"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  plots <- list(
    overall_na_plot = p1,
    na_by_exp_plot = p2,
    sgRNA_duplicates_plot = p3
  )
  
  if (isTRUE(print)) {
    invisible(lapply(plots, print))
  }
  
  return(plots)
}


control_sanity_check <- function(maude_counts_df, cfg, print = TRUE) {
  
  unsorted_bin <- cfg$bins %>%
    dplyr::filter(sorted_or_unsorted == "unsorted") %>%
    dplyr::pull(bin_name) %>%
    base::unname()
  
  sorted_bins <- cfg$bins %>%
    dplyr::filter(sorted_or_unsorted == "sorted") %>%
    dplyr::pull(bin_name)
  
  plot_list <- list()
  
  for (sublib in unique(maude_counts_df$sublib)) {
    
    sublib_data <- maude_counts_df[maude_counts_df$sublib == sublib, ]
    control_data <- sublib_data %>% dplyr::filter(isNontargeting)
    
    if (nrow(control_data) == 0) {
      logger::log_warn("No non-targeting sgRNAs found for sublibrary: {sublib}")
      next
    }
    
    # Devide each bin by the unsorted bin. 
    log2_values <- do.call(
      rbind,
      lapply(
        sorted_bins,
        function(current_bin) {
          data.frame(
            sgRNA = control_data$sgRNA,
            log2_fold_change = log2(
              control_data[[current_bin]] /
                control_data[[unsorted_bin]]
            ),
            comparison = paste(current_bin, "vs", unsorted_bin),
            stringsAsFactors = FALSE
          )
        }
      )
    )
    
    p <- ggplot(
      log2_values,
      aes(
        x = comparison,
        y = log2_fold_change,
        fill = comparison
      )
    ) +
      geom_boxplot(width = 0.5, alpha = 0.7) +
      theme_bw() +
      labs(
        title = paste(
          "Log2 Fold Change Between Sorted Bins and",
          unsorted_bin,
          "for",
          sublib
        ),
        x = "",
        y = "Log2 Fold Change"
      ) +
      theme(
        axis.text.x = element_text(
          angle = 50,
          hjust = 1,
          face = "bold",
          size = 12
        ),
        plot.title = element_text(size = 12),
        legend.position = "none"
      )
    
    safe_sublib <- gsub("[^A-Za-z0-9_.-]+", "_", sublib)
    plot_list[[paste0("control_sanity_", safe_sublib)]] <- p
  }
  
  control_data_all <- maude_counts_df %>% dplyr::filter(isNontargeting)
  
  if (nrow(control_data_all) > 0) {
    
    log2_values_all <- do.call(
      rbind,
      lapply(
        sorted_bins,
        function(current_bin) {
          data.frame(
            sgRNA = control_data_all$sgRNA,
            log2_fold_change = log2(
              control_data_all[[current_bin]] /
                control_data_all[[unsorted_bin]]
            ),
            comparison = paste(current_bin, "vs", unsorted_bin),
            stringsAsFactors = FALSE
          )
        }
      )
    )
    
    p_all <- ggplot(
      log2_values_all,
      aes(
        x = comparison,
        y = log2_fold_change,
        fill = comparison
      )
    ) +
      geom_boxplot(width = 0.5, alpha = 0.7) +
      theme_bw() +
      labs(
        title = paste(
          "Log2 Fold Change Between Sorted Bins and",
          unsorted_bin,
          "for All Data"
        ),
        x = "",
        y = "Log2 Fold Change"
      ) +
      theme(
        axis.text.x = element_text(
          angle = 50,
          hjust = 1,
          face = "bold",
          size = 12
        ),
        plot.title = element_text(size = 12),
        legend.position = "none"
      )
    
    plot_list[["control_sanity_all"]] <- p_all
    
  } else {
    logger::log_warn("No non-targeting sgRNAs found in MAUDE counts.")
  }
  
  if (isTRUE(print)) {
    invisible(lapply(plot_list, print))
  }
  
  return(plot_list)
}


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