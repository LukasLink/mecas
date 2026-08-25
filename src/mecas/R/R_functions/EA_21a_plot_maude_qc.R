# R/R_functions/EA_21a_plot_maude_qc.R


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

