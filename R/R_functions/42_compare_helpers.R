# Helper function to apply FDR filter
apply_FDR_filter <- function(df, FDR_threshold) {
  df <- df %>%
    filter(FDR <= FDR_threshold)
  return(df)
}

# Helper function to calculate number of overlapping entrez
calculate_overlap <- function(reference_df, target_df) {
  reference_entrez <- reference_df$entrez
  target_entrez <- target_df$entrez
  overlap_count <- sum(target_entrez %in% reference_entrez)
  return(overlap_count)
}
# Helper function to calculate the correlation between a refrence_df
# and a target_df
calculate_correlation <- function(reference_df,
                                  target_df,
                                  method = "spearman",
                                  use = "complete.obs",
                                  FDR_threshold = 0.05,
                                  significance_filter = NULL) {
  # keep only needed cols, drop duplicates, sort by FDR, assign ranks
  if (!is.null(significance_filter)) {
    reference_df <- reference_df %>% filter(significanceZ * significance_filter > 0)
    target_df <- target_df %>% filter(significanceZ * significance_filter > 0)
  }
  
  ref <- reference_df[, c("entrez", "FDR")]
  tgt <- target_df[, c("entrez", "FDR")]
  
  ref <- ref[!is.na(ref$entrez) & !is.na(ref$FDR), ]
  tgt <- tgt[!is.na(tgt$entrez) & !is.na(tgt$FDR), ]
  
  ref <- ref %>% rename(ref_FDR = FDR)
  tgt <- tgt %>% rename(tgt_FDR = FDR)
  
  # pair by entrez (inner join)
  paired <- merge(ref, tgt, by = "entrez") %>%
    filter(ref_FDR <= FDR_threshold | tgt_FDR <= FDR_threshold)
  # if not enough paired points, correlation is undefined
  if (nrow(paired) < 2) return(NA_real_)
  
  correlation <- cor(paired$ref_FDR, paired$tgt_FDR, method = method, use = use)
  correlation <- round(correlation, 3)
  return(correlation)
}
# Function to process each scenario
# generates a dataframe with hits, and overlap to reference for each df in df_list

make_compare_df <- function(df_list, name_list, FDR_threshold, significance_filter = NULL, top_N = NULL) {
  
  if (length(df_list) != length(name_list)){
    stop(paste("ERROR: in make_compare_df -> length of df_list:",length(df_list), " must be equal to lenght of name_list:",length(name_list)))
  }
  # First dataframe as the reference (100% by definition)
  reference_df <- df_list[[1]] %>%
    apply_FDR_filter(FDR_threshold)
  if (!is.null(significance_filter)) {
    reference_df <- reference_df %>% filter(significanceZ * significance_filter > 0)
  }
  
  # Initialize empty lists for storing results
  Overlap <- list()
  Overlap_percentage <- list()
  Total_Hits <- list()
  
  # Loop through the dataframes and calculate overlap for each subsample percentage
  for (i in seq_along(df_list)) {
    
    target_df <- df_list[[i]]
    
    # Apply FDR filter to the target dataframe
    target_df <- apply_FDR_filter(target_df, FDR_threshold)
    
    # Apply additional significance filter if provided
    if (!is.null(significance_filter)) {
      target_df <- target_df %>% filter(significanceZ * significance_filter > 0)
    }
    
    # Check overlap with the reference
    overlap_count <- calculate_overlap(reference_df, target_df)
    
    # Calculate percentage overlap
    overlap_percentage_value <- round((overlap_count / nrow(reference_df)) * 100,2)
    
    total_hits_value <- nrow(target_df)
    
    # Append results to the lists
    Overlap <- c(Overlap, overlap_count)
    Overlap_percentage <- c(Overlap_percentage, overlap_percentage_value)
    Total_Hits <- c(Total_Hits, total_hits_value)
  }
  # Combine the lists into a final dataframe
  final_results <- data.frame(
    Name = name_list,
    Overlap = unlist(Overlap),
    Overlap_percentage = unlist(Overlap_percentage),
    Total_Hits = unlist(Total_Hits)
  )
  
  # Return the final dataframe
  return(final_results)
}
make_venn_from_overlap_df <- function(df,
                                      ncol = 3,
                                      main_title = "Venn Diagramm",
                                      comparing = "") {
  stopifnot(all(c("Name", "Overlap", "Total_Hits") %in% colnames(df)))
  if (nrow(df) < 2) stop("df needs at least 2 rows (1 reference + >=1 target).")
  
  if (comparing != ""){
    if (comparing == "subsamples"){
      df <- df %>% 
        mutate(Name = paste0("sub_perc_",Name))
    } else {
      stop("currently only '' (sublib) and 'subsamples' are viable for comparing ")
    }
  }
  
  
  ref_name <- df$Name[1]
  ref_total <- df$Total_Hits[1]
  
  grobs <- list()
  #This is supposed to stop it drawing a venn diagramm with draw.pairwise.venn
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  
  for (i in 2:nrow(df)) {
    tgt_name  <- df$Name[i]
    tgt_total <- df$Total_Hits[i]
    overlap   <- df$Overlap[i]
    
    # sanity (avoid negative/invalid overlaps)
    overlap <- max(0, overlap)
    overlap <- min(overlap, ref_total, tgt_total)
    
    
    g_list <- VennDiagram::draw.pairwise.venn(
      area1 = ref_total,
      area2 = tgt_total,
      cross.area = overlap,
      category = c(ref_name, tgt_name),
      fill = c("lightgreen", "lightblue"),
      alpha = c(0.6, 0.6),
      cat.cex = 1.0,
      cex = 1.3,
      col = c("darkgreen", "darkblue"),
      cat.fontface = "bold",
      cat.col = c("darkgreen", "darkblue"),      # category label colours
      label.col = c("darkgreen", "black", "darkblue"),  # left / overlap / right numbers
      cat.pos  = c(-60, 60),
      cat.dist = c(0.2, 0.2),
      ind = TRUE
    )
    
    # wrap list-of-grobs into a single grobTree so gridExtra can place it
    g <- grid::grobTree(gList = do.call(grid::gList, g_list))
    g <- grid::grobTree(g, vp = grid::viewport(width = 0.5, height = 0.7))
    grobs[[tgt_name]] <- g
  }
  
  title_grob <- grid::textGrob(main_title, gp = grid::gpar(fontsize = 14, fontface = "bold"))
  
  return(
    gridExtra::arrangeGrob(
      title_grob,
      gridExtra::arrangeGrob(grobs = grobs, ncol = ncol),
      ncol = ncol,
      heights = c(0.08, 0.92)
    )
  )
}
# Make the Correlation Heatmap
make_correlation_heatmap <- function(df_list,
                                     name_list,
                                     comparing = "",
                                     plot_title = NULL,
                                     significance_filter = NULL,
                                     FDR_threshold = 0.05) {
  
  if (comparing == "subsamples"){
    name_list = paste0("sub_perc_",name_list)
  }
  
  n <- length(df_list)
  
  # correlation matrix
  cor_mat <- matrix(
    NA_real_,
    nrow = n, ncol = n,
    dimnames = list(name_list, name_list)
  )
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      cor_mat[i, j] <- calculate_correlation(
        df_list[[i]],
        df_list[[j]],
        FDR_threshold = FDR_threshold,
        significance_filter = significance_filter
      )
    }
  }
  
  cor_df <- as.data.frame(as.table(cor_mat))
  colnames(cor_df) <- c("Reference", "Target", "Correlation")
  
  p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Target, y = Reference, fill = Correlation)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(Correlation), "", sprintf("%.3f", Correlation))),
      size = 3
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = paste("Pairwise Spearman Correlation of FDR Ranks for",plot_title),
      x = NULL, y = NULL,
      fill = "Rank corr"
    )
  
  return(p)
}
plot_overlap <- function(df, title, comparing = "") {
  # Ensure the input dataframe contains the necessary columns
  if (!all(c("Name", "Overlap", "Overlap_percentage", "Total_Hits") %in% colnames(df))) {
    stop("The dataframe must contain 'Subsample_percentage', 'Overlap', 'Total_Hits' and 'Overlap_percentage' columns.")
  }
  
  # Create the plot
  if (comparing == ""){
    p <- ggplot(df, aes(x = Name, y = Overlap_percentage)) +
      geom_point(size = 3) +
      theme_bw() +
      labs(
        title = paste("Percentage of correct Hits for", title),
        x = "",
        y = "Overlap Percentage"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
      )  
  }
  if (comparing == "subsamples") {
    
    # Extract trailing number from tags like "A25" / "subX10" / "25"
    df$Number <- as.integer(sub("^\\D*([0-9]+)$", "\\1", df$Name))
    
    # Safety check (optional)
    if (any(is.na(df$Number))) {
      stop("Some 'Name' entries do not end in a number (cannot extract Number).")
    }
    
    # Sort by numeric value so the line goes left -> right correctly
    df <- df[order(df$Number), ]
    
    if (anyDuplicated(df$Number) == 0) {
      # mbers are unique → simple plot
      p <- ggplot(df, aes(x = Number)) +
        geom_line(aes(y = Overlap_percentage), size = 1, group = 1) +
        geom_point(aes(y = Overlap_percentage), size = 3) +
        labs(
          title = paste("Percentage of correct Hits for", title),
          x = "Subsample Percentage",
          y = "Overlap Percentage"
        ) +
        theme_bw()
      
    } else {
      # Numbers NOT unique → plot mean per Number + error bars
      
      # Summarise mean + SE (swap to sd if you prefer)
      tmp <- split(df$Overlap_percentage, df$Number)
      df_sum <- data.frame(
        Number = as.integer(names(tmp)),
        mean   = sapply(tmp, mean),
        sd     = sapply(tmp, sd),
        n      = sapply(tmp, length)
      )
      
      df_sum$se <- df_sum$sd / sqrt(df_sum$n)
      df_sum$se[is.na(df_sum$se)] <- 0  # happens when n=1 (sd=NA)
      
      df_sum <- df_sum[order(df_sum$Number), ]
      
      p <- ggplot(df_sum, aes(x = Number)) +
        geom_line(aes(y = mean), size = 1, group = 1) +
        geom_point(aes(y = mean), size = 3) +
        geom_errorbar(
          aes(ymin = mean - se, ymax = mean + se),
          width = 0.2
        ) +
        labs(
          title = paste("Percentage of correct Hits with SE for", title),
          x = "Subsample Percentage",
          y = "Overlap Percentage"
        ) +
        theme_bw()
      
      # Optional: show the individual points too (nice for debugging)
      # p <- p + geom_point(data = df, aes(x = Number, y = Overlap_percentage), alpha = 0.4)
    }
  }
  
  return(p)
}

# Function for a line plot comparing total hits
# Currently depricated and replaced by venn diagramms
plot_total_hits <- function(df, title, comparing = "") {
  # Ensure the input dataframe contains the necessary columns
  if (!all(c("Name", "Overlap", "Overlap_percentage", "Total_Hits") %in% colnames(df))) {
    stop("The dataframe must contain 'Name', 'Overlap', 'Total_Hits' and 'Overlap_percentage' columns.")
  }
  if (comparing == ""){
    x_lab_text = ""
  }
  if (comparing == "subsamples"){
    x_lab_text = "Subsample Percentage"
  }
  # Create the plot
  p <- ggplot(df, aes(x = Name)) +
    geom_line(aes(y = Total_Hits), size = 1) + 
    geom_point(aes(y = Total_Hits), size = 3) +
    labs(
      title = paste("Number of all Hits (not necessarily correct ones) for",title),
      x = x_lab_text,
      y = "Total Hits"
    ) +
    theme_bw()
  
  # Return the plot
  return(p)
}

compare_df_list <- function(df_list,
                            name_list,
                            FDR_threshold = 0.05,
                            top_N = NULL,
                            comparing = "",
                            correlation_heatmap = TRUE,
                            overlap = TRUE,
                            venn_diagram = TRUE) {
  
  
  all_hits_df <- make_compare_df(df_list, name_list, FDR_threshold)
  print(all_hits_df)
  sig_up_hits_df <- make_compare_df(df_list, name_list, FDR_threshold, significance_filter = 1)
  sig_down_hits_df <- make_compare_df(df_list, name_list, FDR_threshold, significance_filter = -1)
  results_list <- list(all_hits_df, sig_up_hits_df, sig_down_hits_df)
  names(results_list) <- c("All_Hits", "Up_Hits" ,"Down_Hits")
  signif_map <- list(
    All_Hits  = NULL,
    Up_Hits   = 1,
    Down_Hits = -1
  )
  
  # Loop over all dataframes in the results list and plot them
  for (scenario_name in names(results_list)) {
    signif_filter_for_scenario <- signif_map[[scenario_name]]
    if (correlation_heatmap == TRUE){
      q <- make_correlation_heatmap(df_list,
                                    name_list,
                                    comparing = comparing,
                                    plot_title = scenario_name,
                                    significance_filter = signif_filter_for_scenario,
                                    FDR_threshold = FDR_threshold)
      print(q)      
    }
    if (overlap == TRUE){
      p <- plot_overlap(results_list[[scenario_name]],
                        scenario_name,
                        comparing = comparing)
      print(p)       
    }
    if (venn_diagram == TRUE){
      r <- make_venn_from_overlap_df(results_list[[scenario_name]],
                                     ncol = 1,
                                     comparing = comparing,
                                     main_title = scenario_name)
      grid::grid.newpage()
      grid::grid.draw(r)      
    }
  }
  return(results_list)
}