# R/R_functions/EA_12_plot_violin_by_sublib_sample.R


plot_violin_by_sublib_sample <- function(count_df_long,
                                         cfg,
                                         norm_method = NULL) {
  df <- count_df_long
  
  if (!is.null(norm_method)) {
    df <- normalize_count_df_long(count_df_long, norm_method)
    norm_title <- paste0("(normalization: ", norm_method, ")")
  } else {
    norm_title <- "(no normalization)"
  }
  
  count_type <- cfg$counting$data_type
  
  if (!(count_type %in% c("umis", "reads"))) {
    stop(
      "Invalid data_type: ", count_type,
      ". data_type must be 'umis' or 'reads'.",
      call. = FALSE
    )
  }
  
  if (count_type == "umis") {
    title_prefix <- "UMI counts for"
  }
  
  if (count_type == "reads") {
    title_prefix <- "Read counts for"
  }
  
  df <- df %>%
    dplyr::mutate(
      group_type = ifelse(group_category == "targeting", "targeting", "non-targeting")
    )
  
  unique_groups <- df %>%
    dplyr::distinct(sublib, sample)
  
  plot_list <- list()
  
  for (i in 1:nrow(unique_groups)) {
    sublib_val <- unique_groups$sublib[i]
    sample_val <- unique_groups$sample[i]
    
    df_subset <- df %>%
      dplyr::filter(sublib == sublib_val, sample == sample_val)
    
    df_subset <- df_subset %>%
      dplyr::rename(Type = group_type) %>% 
      dplyr::group_by(Type, bin_name) %>%
      dplyr::filter(dplyr::n() > 0) %>%
      dplyr::ungroup()
    
    p <- ggplot2::ggplot(
      df_subset,
      ggplot2::aes(x = bin_name, y = count, fill = Type)
    ) +
      ggplot2::geom_violin(trim = FALSE, scale = "width") +
      ggplot2::labs(
        title = paste0(title_prefix, " ", sublib_val, ", ", sample_val, " ", norm_title),
        x = "",
        y = "Count"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    plot_list[[paste(sublib_val, sample_val, sep = "_")]] <- p
  }
  
  return(plot_list)
}


