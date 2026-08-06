# R/R_functions/EA_21b_control_sanity_check.R
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

