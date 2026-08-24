# R/R_functions/DA_13_normalize_count_df_long.R

normalize_count_df_long <- function(count_df_long,
                                    umis_as_sublibs = FALSE,
                                    norm_method = "control_median",
                                    return_info = FALSE) {
  allowed_norm_methods <- c("control_median")
  
  if (!(norm_method %in% allowed_norm_methods)) {
    cat("ERROR: ", norm_method, "is not an allowed normalization method.\n")
    cat("Implemented normalization methods: ", allowed_norm_methods, "\n")
    stop()
  }
  
  pseudocount_added <- FALSE
  
  if (norm_method == "control_median") {
    log_info("Performing Normalization for each input file (pair) using the median of the control sgRNAs.")
    
    # Define Grouping Columns for our two cases: 
    # 1. Sublibs define seperately sequenced samples
    # 2. Sublibs define a UMI
    if (isTRUE(umis_as_sublibs)) {
      norm_group_cols <- c("bin_name", "sample")
    } else {
      norm_group_cols <- c("bin_name", "sublib", "sample")
    }
    
    norm_fac <- count_df_long %>%
      dplyr::filter(group_category %in% c(
        "targeting_control",
        "non_targeting_control",
        "kept_control"
      )) %>%
      dplyr::group_by(
        dplyr::across(dplyr::all_of(norm_group_cols))
      ) %>%
      dplyr::summarise(norm_factor = median(count), .groups = "drop")
    
    med_count <- median(count_df_long$count)
    
    count_df_long_continue <- count_df_long
    
    if (med_count == 0 | any(norm_fac$norm_factor == 0)) {
      cat("WARNING: at least one normalization factor is 0\n")
      cat("Median of all counts:", med_count, "\n")
      
      zero_nf <- norm_fac %>%
        dplyr::filter(norm_factor == 0)
      
      if (!(nrow(zero_nf) == 0)) {
        cat("Groups with norm_factor == 0:\n")
        print(norm_fac)
        cat("Adding global pseudocount (+1).\n")
        
        pseudocount_added <- TRUE
        
        count_df_long_plus_one <- count_df_long %>%
          dplyr::mutate(count = count + 1)
        
        norm_fac <- count_df_long_plus_one %>%
          dplyr::filter(group_category %in% c(
            "targeting_control",
            "non_targeting_control",
            "kept_control"
          )) %>%
          dplyr::group_by(
            dplyr::across(dplyr::all_of(norm_group_cols))
          ) %>%
          dplyr::summarise(norm_factor = median(count), .groups = "drop")
        
        med_count <- median(count_df_long_plus_one$count)
        
        count_df_long_continue <- count_df_long_plus_one
      }
      
      cat("--------------------------------------------\n")
    }
    
    return_df <- count_df_long_continue %>%
      dplyr::inner_join(norm_fac, by = norm_group_cols) %>%
      dplyr::mutate(norm_count = (count * med_count) / norm_factor) %>%
      dplyr::mutate(count = round(norm_count, 2)) %>%
      dplyr::select(-c(norm_factor, norm_count))
  }
  
  return_df <- return_df %>%
    dplyr::mutate(count = round(count))
  
  if (isTRUE(return_info)) {
    return(list(
      count_df_long = return_df,
      pseudocount_already_added = pseudocount_added
    ))
  }
  
  return(return_df)
}