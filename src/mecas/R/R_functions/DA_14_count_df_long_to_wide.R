# R/R_functions/DA_14_count_df_long_to_wide.R


count_df_long_to_wide <- function(count_df_long,
                                  bins,
                                  print = TRUE,
                                  drop_0s = TRUE,
                                  recover_input = FALSE,
                                  for_sum = FALSE,
                                  pseudocount_already_added = FALSE,
                                  strict_mode = FALSE,
                                  umis_as_sublibs = FALSE, 
                                  ctrl_selection = c("targeting_control", "non_targeting_control")) {
  
  required_count_columns <- c("sgRNA", "count", "group_category", "bin_name", "sublib", "sample", "exp")
  missing_count_columns <- setdiff(required_count_columns, colnames(count_df_long))
  
  if (length(missing_count_columns) > 0) {
    stop_log(
      "`count_df_long` is missing required columns:\n",
      paste0("  - ", missing_count_columns, collapse = "\n")
    )
  }
  
  unsorted_bin <- bins %>%
    dplyr::filter(sorted_or_unsorted == "unsorted") %>%
    dplyr::pull(bin_name)
  
  sorted_bins <- bins %>%
    dplyr::filter(sorted_or_unsorted == "sorted") %>%
    dplyr::arrange(bin_fraction_min, bin_fraction_max) %>%
    dplyr::pull(bin_name)
  
  all_bins <- c(unsorted_bin, sorted_bins)
  
  log_wide_diagnostics <- function(df, stage) {
    if (!isTRUE(print)) {
      return(invisible(NULL))
    }
    
    logger::log_info("{stage} dimension: {paste(dim(df), collapse = ' x ')}")
    logger::log_info("{stage} total NAs: {sum(is.na(df))}")
    
    for (current_bin in all_bins) {
      logger::log_info("{stage} NAs in `{current_bin}`: {sum(is.na(df[[current_bin]]))}")
    }
    
    invisible(NULL)
  }
  
  if (isTRUE(print)) {
    logger::log_info("NAs in count_df_long: {sum(is.na(count_df_long))}")
  }
  
  maude_counts_df <- count_df_long %>%
    dplyr::mutate(isNontargeting = group_category %in% ctrl_selection) %>%
    dplyr::select(-group_category) %>%
    tidyr::pivot_wider(names_from = bin_name, values_from = count)
  
  # If a guide is present in at least one sorted bin, missing values in the
  # other sorted bins represent zero observed counts.
  has_any_sorted_count <- rowSums(!is.na(maude_counts_df[, sorted_bins, drop = FALSE])) > 0
  missing_count_value <- if (isTRUE(pseudocount_already_added)) 1 else 0
  
  maude_counts_df <- maude_counts_df %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(sorted_bins),
        ~ ifelse(is.na(.x) & has_any_sorted_count, missing_count_value, .x)
      )
    )
  
  log_wide_diagnostics(maude_counts_df, "Initial wide table")
  
  #-----------------------------------------------------------------------------
  # Optional: recover input (by estimating the mean of other bins for each sgRNA)
  #-----------------------------------------------------------------------------
  if (isTRUE(recover_input)) {
    if (isTRUE(umis_as_sublibs)){
      log_warn("recover_input option is currently not supported for --umis-as-sublibs and will be skipped")
    } else {
      input_means <- maude_counts_df %>%
        dplyr::group_by(sgRNA, sublib) %>%
        dplyr::summarise(mean_input = mean(.data[[unsorted_bin]], na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(mean_input = ifelse(is.nan(mean_input), NA_real_, mean_input))
      
      maude_counts_df <- maude_counts_df %>%
        dplyr::left_join(input_means, by = c("sgRNA", "sublib"))
      
      recover_input_rows <- is.na(maude_counts_df[[unsorted_bin]]) & has_any_sorted_count
      maude_counts_df[[unsorted_bin]][recover_input_rows] <- maude_counts_df$mean_input[recover_input_rows]
      maude_counts_df$mean_input <- NULL
    }
  } else if (!isTRUE(for_sum)) {
    maude_counts_df <- maude_counts_df %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(all_bins), ~ !is.na(.x)))
  }
  #-----------------------------------------------------------------------------
  
  if (isTRUE(drop_0s)) {
    maude_counts_df <- maude_counts_df %>%
      dplyr::filter(!dplyr::if_all(dplyr::all_of(all_bins), ~ .x == 0))
  }
  
  if (!isTRUE(for_sum) && !isTRUE(pseudocount_already_added)) {
    maude_counts_df <- maude_counts_df %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(all_bins), ~ .x + 1))
  }
  
  #-----------------------------------------------------------------------------
  # Optional: strict mode
  #-----------------------------------------------------------------------------
  if (isTRUE(strict_mode)) {
    logger::log_info("Strict mode enabled.")
    logger::log_info("Rows/Guides before strict mode: {nrow(maude_counts_df)}")
    
    maude_counts_df <- maude_counts_df %>%
      dplyr::filter(
        dplyr::if_all(
          dplyr::all_of(all_bins),
          ~ !is.na(.x) & .x > 1
        )
      )
    logger::log_info("Rows/Guides after strict mode: {nrow(maude_counts_df)}")
  }
  #-----------------------------------------------------------------------------
  
  # MAUDE expects integer counts.
  maude_counts_df <- maude_counts_df %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(all_bins), round)) %>%
    as.data.frame()
  
  log_wide_diagnostics(maude_counts_df, "Final wide table")
  
  return(maude_counts_df)
}