# These functions are for creating violin plots for initial data inspection.
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

get_grouped_summary_wide <- function(count_df_long,
                                     stat = c("median", "mean")) {
  
  stat <- match.arg(stat)
  
  summary_fun <- switch(
    stat,
    median = function(x) median(x, na.rm = TRUE),
    mean = function(x) round(mean(x, na.rm = TRUE), 2)
  )
  
  count_df_long %>%
    dplyr::mutate(
      group_type = ifelse( group_category == "targeting", "targeting", "non-targeting")
    ) %>%
    dplyr::group_by(sublib, sample, group_type, bin_name) %>%
    dplyr::summarise(value = summary_fun(count), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = bin_name, values_from = value) %>%
    dplyr::mutate(group_id = paste(sublib, sample, group_type, sep = "_")) %>%
    dplyr::select(group_id, dplyr::everything(), -sublib, -sample, -group_type)
}

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

generate_all_violin_plots_and_summaries <- function(df,
                                                    cfg,
                                                    norm_method = NULL,
                                                    targeting = TRUE,
                                                    non_targeting = FALSE,
                                                    summary_df = FALSE,
                                                    y_limit = 60) {
  if (isTRUE(targeting)) {
    targeting_results <- plot_violin_by_group_category_split_by_sublib(
      df,
      cfg = cfg,
      include_targeting = TRUE,
      norm_method = norm_method,
      y_limit = y_limit,
      viol_col = "lightgreen"
    )
    
    cat("Targeting plots:\n")
    
    for (sublib in names(targeting_results$plots)) {
      print(targeting_results$plots[[sublib]])
    }
    
    if (summary_df == TRUE) {
      cat("\nTargeting count summary:\n")
      print(targeting_results$count_summary)
    }
  }
  
  if (isTRUE(non_targeting)) {
    non_targeting_results <- plot_violin_by_group_category_split_by_sublib(
      df,
      cfg = cfg,
      include_targeting = FALSE,
      norm_method = norm_method,
      y_limit = y_limit,
      viol_col = "lightblue"
    )
    
    cat("\nNon-targeting plots:\n")
    
    for (sublib in names(non_targeting_results$plots)) {
      print(non_targeting_results$plots[[sublib]])
    }
    
    if (summary_df == TRUE) {
      cat("\nNon-targeting count summary:\n")
      print(non_targeting_results$count_summary)
    }
  }
}