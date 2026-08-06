# R/R_functions/EA_14_generate_all_violin_plots_and_summaries.R

# I think this function is deprecated now. 

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