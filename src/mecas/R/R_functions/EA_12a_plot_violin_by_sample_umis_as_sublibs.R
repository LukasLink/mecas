# R/R_functions/EA_12a_plot_violin_by_sample_umis_as_sublibs.R


plot_violin_by_sample_umis_as_sublibs <- function(count_df_long,
                                         cfg,
                                         norm_method = NULL) {
  df <- count_df_long
  
  if (!is.null(norm_method)) {
    df <- normalize_count_df_long(count_df_long, norm_method, umis_as_sublibs = TRUE)
    norm_title <- paste0("(normalization: ", norm_method, ")")
  } else {
    norm_title <- "(no normalization)"
  }

  df <- df %>%
    dplyr::mutate(
      Type = ifelse(group_category == "targeting", "targeting", "non-targeting")
    )
  
  plot_list <- list()
  
  for (sample_val in unique(df$sample)) {
    
    df_subset <- df %>%
      dplyr::filter(sample == sample_val)
    
    p <- ggplot2::ggplot(
      df_subset,
      ggplot2::aes(x = bin_name, y = count, fill = Type)
    ) +
      ggplot2::geom_violin(trim = FALSE, scale = "width") +
      ggplot2::labs(
        title = paste0("Reads per UMI for ",sample_val," ",norm_title),
        x = "",
        y = "Reads per UMI"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45,hjust = 1)
      )
    
    plot_list[[sample_val]] <- p
  }
  
  return(plot_list)
}


