# R/R_functions/EA_41_scatter_plot_for_one_group.R

scatter_plot_for_one_group <- function(df) {
  
  make_safe_name <- function(x) {
    x <- ifelse(is.na(x) | x == "", "NA", x)
    gsub("[^A-Za-z0-9_.-]+", "-", x)
  }
  
  sample_name <- unique(df$sample)[1]
  bin_name    <- unique(df$bin_name)[1]
  sublib      <- unique(df$sublib)[1]
  
  title_txt <- paste(
    "Read count vs UMI count",
    paste0("sample = ", sample_name),
    paste0("bin_name = ", bin_name),
    paste0("sublib = ", sublib),
    sep = "\n"
  )
  
  cor_val <- suppressWarnings(
    cor(df$read_count, df$umi_count, method = "spearman", use = "complete.obs")
  )
  
  p <- ggplot(df, aes(x = read_count, y = umi_count)) +
    geom_point(alpha = 0.35, size = 1) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = title_txt,
      x = "Read count +1 (log10 scale)",
      y = "UMI count +1 (log10 scale)"
    ) +
    theme_bw()
  
  if (is.finite(cor_val)) {
    p <- p +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = paste0("Spearman rho = ", round(cor_val, 3)),
        hjust = 1.1, vjust = 1.5, size = 4
      )
  }
  
  list(
    plot = p,
    filename = paste0(
      "umi_read_scatter_",
      make_safe_name(sample_name), "_",
      make_safe_name(bin_name), "_",
      make_safe_name(sublib),
      ".png"
    )
  )
}