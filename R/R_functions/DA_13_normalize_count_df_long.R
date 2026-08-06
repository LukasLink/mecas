# R/R_functions/DA_13_normalize_count_df_long.R

normalize_count_df_long <- function(count_df_long,
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
    
    norm_fac <- count_df_long %>%
      dplyr::filter(group_category %in% c(
        "targeting_control",
        "non_targeting_control",
        "kept_control"
      )) %>%
      dplyr::group_by(bin_name, sublib, sample) %>%
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
          dplyr::group_by(bin_name, sublib, sample) %>%
          dplyr::summarise(norm_factor = median(count), .groups = "drop")
        
        med_count <- median(count_df_long_plus_one$count)
        
        count_df_long_continue <- count_df_long_plus_one
      }
      
      cat("--------------------------------------------\n")
    }
    
    return_df <- count_df_long_continue %>%
      dplyr::inner_join(norm_fac, by = c("bin_name", "sublib", "sample")) %>%
      dplyr::mutate(norm_count = (count * med_count) / norm_factor) %>%
      dplyr::mutate(count = round(norm_count, 2)) %>%
      dplyr::select(-c(norm_factor, norm_count))
  }
  
  return_df <- return_df %>%
    dplyr::mutate(count = round(count))
  
  if (isTRUE(return_info)) {
    return(list(
      count_df_long = return_df,
      pseudocount_added = pseudocount_added
    ))
  }
  
  return(return_df)
}