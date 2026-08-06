# R/R_functions/EA_12_plot_violin_by_group_category_split_by_sublib.R


plot_violin_by_group_category_split_by_sublib <- function(df,
                                                          cfg,
                                                          include_targeting = TRUE,
                                                          norm_method = NULL,
                                                          y_limit = 60,
                                                          box_col = "white",
                                                          viol_col = "#E0E0E0") {
  
  if (!is.null(norm_method)) {
    df <- normalize_count_df_long(df, norm_method)
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
  
  filtered_df <- df %>%
    dplyr::filter(
      if (include_targeting) group_category == "targeting"
      else group_category != "targeting"
    ) %>%
    dplyr::mutate(group = interaction(bin_name, exp, drop = TRUE))
  
  count_summary <- filtered_df %>%
    dplyr::group_by(bin_name, exp) %>%
    dplyr::summarise(
      total_count = sum(count),
      mean_count = round(mean(count), 1),
      sd_count = round(stats::sd(count), 1),
      .groups = "drop"
    ) %>%
    dplyr::mutate(group = interaction(bin_name, exp, drop = TRUE))
  
  plots <- list()
  
  for (sublib_name in unique(filtered_df$sublib)) {
    sub_df <- filtered_df %>%
      dplyr::filter(sublib == sublib_name)
    
    x_levels <- levels(interaction(sub_df$bin_name, sub_df$exp, drop = TRUE))
    
    label_df <- count_summary %>%
      dplyr::filter(
        interaction(bin_name, exp, drop = TRUE) %in%
          unique(interaction(sub_df$bin_name, sub_df$exp, drop = TRUE))
      ) %>%
      dplyr::mutate(
        group_fac = factor(interaction(bin_name, exp, drop = TRUE), levels = x_levels),
        x_center = as.numeric(group_fac)
      )
    
    p <- ggplot2::ggplot(
      sub_df,
      ggplot2::aes(x = interaction(bin_name, exp), y = count)
    ) +
      ggplot2::geom_violin(
        fill = viol_col,
        color = NA,
        scale = "width",
        trim = FALSE
      ) +
      ggplot2::geom_boxplot(
        width = 0.1,
        outlier.size = 0.2,
        fill = box_col
      ) +
      ggplot2::geom_rect(
        data = label_df,
        ggplot2::aes(
          xmin = x_center - 0.3,
          xmax = x_center + 0.3,
          ymin = y_limit - y_limit / 8 - 0.1,
          ymax = y_limit + 0.1
        ),
        inherit.aes = FALSE,
        fill = "white",
        color = "black",
        linewidth = 0.15
      ) +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = group, y = y_limit - y_limit / 26, label = total_count),
        size = 2.5,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = group, y = y_limit - y_limit / 14, label = paste0("mean: ", mean_count)),
        size = 2.3,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = group, y = y_limit - y_limit / 10, label = paste0("sd: ", sd_count)),
        size = 2.3,
        inherit.aes = FALSE
      ) +
      ggplot2::labs(
        title = paste(
          title_prefix,
          ifelse(include_targeting, "targeting sgRNA", "non-targeting sgRNA"),
          "–",
          sublib_name,
          norm_title
        ),
        x = "",
        y = "Count"
      ) +
      ggplot2::coord_cartesian(ylim = c(0, y_limit)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    plots[[sublib_name]] <- p
  }
  
  return(list(
    plots = plots,
    count_summary = count_summary
  ))
}
