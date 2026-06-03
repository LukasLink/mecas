plot_maude_qc <- function(maude_counts_df, print = TRUE, add_0s = FALSE) {
  
  if (isTRUE(add_0s)) {
    maude_counts_df <- maude_counts_df %>%
      mutate(
        upper = ifelse(is.na(upper), 0, upper),
        lower = ifelse(is.na(lower), 0, lower)
      )
  }
  
  na_summary_overall <- maude_counts_df %>%
    pivot_longer(cols = c(input, lower, upper), names_to = "condition", values_to = "value") %>%
    mutate(is_na = is.na(value)) %>%
    group_by(condition) %>%
    summarise(na_count = sum(is_na), .groups = "drop")
  
  p1 <- ggplot(na_summary_overall, aes(x = condition, y = na_count, fill = condition)) +
    geom_col() +
    labs(title = "Missing Entries per Bin", y = "Number of NAs", x = "Condition") +
    theme_bw()
  
  na_summary_by_exp <- maude_counts_df %>%
    pivot_longer(cols = c(input, lower, upper), names_to = "condition", values_to = "value") %>%
    mutate(is_na = is.na(value)) %>%
    group_by(exp, condition) %>%
    summarise(na_count = sum(is_na), .groups = "drop")
  
  p2 <- ggplot(na_summary_by_exp, aes(x = exp, y = na_count, fill = condition)) +
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


control_sanity_check <- function(maude_counts_df, print = TRUE) {
  
  plot_list <- list()
  
  for (sublib in unique(maude_counts_df$sublib)) {
    
    sublib_data <- maude_counts_df[maude_counts_df$sublib == sublib, ]
    control_data <- sublib_data[grepl("^CONTROL_", sublib_data$sgRNA), ]
    
    if (nrow(control_data) == 0) {
      logger::log_warn("No CONTROL_ sgRNAs found for sublibrary: {sublib}")
      next
    }
    
    control_data$log2_upper_input <- log2(control_data$upper / control_data$input)
    control_data$log2_lower_input <- log2(control_data$lower / control_data$input)
    
    log2_values <- data.frame(
      sgRNA = rep(control_data$sgRNA, 2),
      log2_fold_change = c(control_data$log2_upper_input, control_data$log2_lower_input),
      comparison = rep(c("upper vs input", "lower vs input"), each = nrow(control_data))
    )
    
    p <- ggplot(log2_values, aes(x = comparison, y = log2_fold_change, fill = comparison)) +
      geom_boxplot(width = 0.5, alpha = 0.7) +
      scale_fill_manual(values = c("upper vs input" = "lightblue", "lower vs input" = "lightgreen")) +
      theme_bw() +
      labs(
        title = paste("Log2 Fold Change Between Upper/Lower and Input for", sublib),
        x = "",
        y = "Log2 Fold Change"
      ) +
      theme(
        axis.text.x = element_text(angle = 50, hjust = 1, face = "bold", size = 12),
        plot.title = element_text(size = 12),
        legend.position = "none"
      )
    
    safe_sublib <- gsub("[^A-Za-z0-9_.-]+", "_", sublib)
    plot_list[[paste0("control_sanity_", safe_sublib)]] <- p
  }
  
  control_data_all <- maude_counts_df[grepl("^CONTROL_", maude_counts_df$sgRNA), ]
  
  if (nrow(control_data_all) > 0) {
    
    control_data_all$log2_upper_input <- log2(control_data_all$upper / control_data_all$input)
    control_data_all$log2_lower_input <- log2(control_data_all$lower / control_data_all$input)
    
    log2_values_all <- data.frame(
      sgRNA = rep(control_data_all$sgRNA, 2),
      log2_fold_change = c(control_data_all$log2_upper_input, control_data_all$log2_lower_input),
      comparison = rep(c("upper vs input", "lower vs input"), each = nrow(control_data_all))
    )
    
    p_all <- ggplot(log2_values_all, aes(x = comparison, y = log2_fold_change, fill = comparison)) +
      geom_boxplot(width = 0.5, alpha = 0.7) +
      scale_fill_manual(values = c("upper vs input" = "lightblue", "lower vs input" = "lightgreen")) +
      theme_bw() +
      labs(
        title = "Log2 Fold Change Between Upper/Lower and Input for All Data",
        x = "",
        y = "Log2 Fold Change"
      ) +
      theme(
        axis.text.x = element_text(angle = 50, hjust = 1, face = "bold", size = 12),
        plot.title = element_text(size = 12),
        legend.position = "none"
      )
    
    plot_list[["control_sanity_all"]] <- p_all
    
  } else {
    logger::log_warn("No CONTROL_ sgRNAs found in MAUDE counts.")
  }
  
  if (isTRUE(print)) {
    invisible(lapply(plot_list, print))
  }
  
  return(plot_list)
}


create_maude_qc_plot_objects <- function(count_df_long,
                                         maude_counts_df,
                                         opt,
                                         input_recovery = FALSE) {
  
  if (identical(opt$data_type, "umis")) {
    data_name <- "UMIs"
  } else if (identical(opt$data_type, "reads")) {
    data_name <- "reads"
  } else {
    stop_log("`data_type` must be either 'umis' or 'reads'.")
  }
  
  if (isTRUE(input_recovery)) {
    title_suffix_1 <- "(without input recovery)"
    title_suffix_2 <- "(with input recovery)"
  } else {
    title_suffix_1 <- "(before MAUDE preparation)"
    title_suffix_2 <- "(ready for MAUDE)"
  }
  
  sample_sum_df <- count_df_long %>%
    group_by(condition, exp) %>%
    summarise(total_count = sum(count), .groups = "drop") %>%
    mutate(Sample = interaction(condition, exp, drop = TRUE),
           Bin = condition)
  
  p_sum_per_sample <- ggplot(sample_sum_df, aes(x = Sample, y = total_count, fill = Bin)) +
    geom_col() +
    labs(
      title = paste("Sum of", data_name, "per sample"),
      x = "Sample",
      y = paste("Sum of", data_name)
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.position = "none"
    )
  
  sample_n_df <- count_df_long %>%
    group_by(condition, exp) %>%
    summarise(num_sgRNAs = n(), .groups = "drop") %>%
    mutate(Sample = interaction(condition, exp, drop = TRUE),
           Bin = condition)
  
  p_n_per_sample <- ggplot(sample_n_df, aes(x = Sample, y = num_sgRNAs, fill = Bin)) +
    geom_col() +
    labs(
      title = paste("Number of sgRNA entries per sample", title_suffix_1),
      x = "Sample",
      y = "Number of sgRNA entries"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.position = "none"
    )
  
  bin_n_before_df <- count_df_long %>%
    group_by(condition) %>%
    summarise(num_sgRNAs = n(), .groups = "drop") %>%
    mutate(Bin = condition)
  
  p_n_per_bin_before <- ggplot(bin_n_before_df, aes(x = Bin, y = num_sgRNAs, fill = Bin)) +
    geom_col() +
    labs(
      title = paste("Number of sgRNA entries per bin", title_suffix_1),
      x = "Bin",
      y = "sgRNA entries"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.position = "none"
    )
  
  bin_n_after_df <- maude_counts_df %>%
    summarise(
      input = sum(!is.na(input)),
      upper = sum(!is.na(upper)),
      lower = sum(!is.na(lower))
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Bin",
      values_to = "count_non_na"
    )
  
  p_n_per_bin_after <- ggplot(bin_n_after_df, aes(x = Bin, y = count_non_na, fill = Bin)) +
    geom_col() +
    labs(
      title = paste("Number of sgRNA entries per bin", title_suffix_2),
      x = "Bin",
      y = "Number of sgRNA entries"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.position = "none"
    )
  
  bin_sum_after_df <- maude_counts_df %>%
    summarise(
      input = sum(input, na.rm = TRUE),
      upper = sum(upper, na.rm = TRUE),
      lower = sum(lower, na.rm = TRUE)
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Bin",
      values_to = "total_count"
    )
  
  p_sum_per_bin_after <- ggplot(bin_sum_after_df, aes(x = Bin, y = total_count, fill = Bin)) +
    geom_col() +
    labs(
      title = paste("Sum of", data_name, "per bin"),
      x = "Bin",
      y = "Total counts"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
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
    maude_counts_df = maude_counts_df,
    print = FALSE
  )
  
  maude_qc_plots <- plot_maude_qc(
    maude_counts_df = maude_counts_df,
    print = FALSE
  )
  
  return(c(
    manual_plots,
    control_plots,
    maude_qc_plots
  ))
}