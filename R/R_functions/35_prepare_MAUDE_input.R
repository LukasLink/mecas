#-----------------------------------------------------------------------------
# prepare data for MAUDE
#-----------------------------------------------------------------------------
run_prepare_data_for_MAUDE <- function(count_df_long,
                                       opt){
  logger::log_info("Preparing Data for MAUDE...")
  if (opt$subsample_controls == TRUE){
    # count_df_long_old <- subsample_controls_func_old(count_df_long, merged_sgRNA_df)
    count_df_long <- subsample_controls_func(count_df_long, merged_sgRNA_df)
  }
  
  if (!(norm_method %in% c("","control_median"))){
    stop("Error: 'norm_method' must be one of '', 'control_median'. The script will now stop.")
  }
  if (norm_method == "control_median"){
    count_df_long <- normalize_count_df_long(count_df_long,
                                             norm_method = norm_method)
  }
  
  
  
  maude_counts_df <- count_df_long_to_wide(count_df_long = count_df_long,
                                           print = FALSE,
                                           drop_0s = drop_0s,
                                           recover_input = recover_input)
  if (!(method %in% c("","rep","sum",'rep_sample','rep_sublib'))){
    stop("Error: 'method' must be one of '', 'rep', 'rep_sample', 'rep_sublib', or 'sum'. The script will now stop.")
  }
  
  if (method == ""){
    maude_counts_df <- maude_counts_df %>% 
      mutate(exp = "rep1")
  }
  if (method == "rep"){
    
  }
  if (method == "rep_sample"){
    maude_counts_df <- maude_counts_df %>%
      mutate(exp = sample)
  }
  if (method == "rep_sublib"){
    maude_counts_df <- maude_counts_df %>%
      mutate(exp = sublib)
  }
  if (method == "sum"){
    stop("Method 'sum' is deprecated, use the option combine_for_guide_stats instead")
    maude_counts_df <- count_df_long_to_wide(count_df_long = count_df_long,
                                             print = FALSE,
                                             drop_0s = drop_0s,
                                             recover_input = TRUE,
                                             for_sum = TRUE)
    # Group by sgRNA and summarize the required columns
    maude_counts_df <- maude_counts_df %>%
      group_by(sgRNA, sublib) %>%
      summarize(
        input = pmax(sum(input, na.rm = TRUE), 0),  # Sum and ensure minimum is 0
        upper = pmax(sum(upper, na.rm = TRUE), 0),  # Sum and ensure minimum is 0
        lower = pmax(sum(lower, na.rm = TRUE), 0),  # Sum and ensure minimum is 0
        isNontargeting = dplyr::first(isNontargeting),  # Take the first value of isNontargeting (same for all in the group)
        .groups = 'drop'  # Drop the group structure after summarizing
      ) %>%
      mutate(
        exp = "rep1",
        input = input + 1,
        upper = upper + 1,
        lower = lower + 1
      )
  }
  
  if (strict_mode){
    if (pseudocount_added){
      umi_threshold <- 2
    } else {
      umi_threshold <- 1
    } 
    cat("Strict Mode enabled\n")
    cat("Rows before strict mode: \t", nrow(maude_counts_df),"\n")
    maude_counts_df <- maude_counts_df %>%
      filter(if_any(c(input, upper, lower), ~ . <= umi_threshold))
    cat("Rows after strict mode: \t", nrow(maude_counts_df),"\n")
  }
  if (length(include_controls_list) > 0) {
    for (control_gene in include_controls_list) {
      maude_counts_df$isNontargeting[grepl(control_gene, maude_counts_df$sgRNA)] <- FALSE
    }
  }
  if (exists("use_only_these_controls_list")) {
    if (length(use_only_these_controls_list) > 0){
      maude_counts_df$isNontargeting[ !(maude_counts_df$sgRNA %in% use_only_these_controls_list) ] <- FALSE
    }
  }
  return(maude_counts_df)
}

#-----------------------------------------------------------------------------
# Pre MAUDE plots
#-----------------------------------------------------------------------------
