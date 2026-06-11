# 05_prepare_data_for_MAUDE

subsample_controls_func_old <- function(count_df_long,
                                        cfg,
                                        merged_sgRNA_df = cfg$merged_sgRNA_df,
                                        percentage_kept_controls = 0.1) {
  
  percentage_used_controls <- 1 - percentage_kept_controls
  
  count_df_long_merged <- merge(
    count_df_long,
    merged_sgRNA_df[, c("sgrna_id", "entrez")],
    by.x = "sgRNA",
    by.y = "sgrna_id",
    all.x = TRUE
  )
  
  count_df_long_merged <- count_df_long_merged %>%
    dplyr::filter(group_category != "targeting") %>%
    dplyr::distinct(sgRNA) %>% 
    dplyr::mutate(
      unique_names = dplyr::case_when(
        grepl("CONTROL_C_[A-Za-z0-9]+_", sgRNA) &
          !grepl("CONTROL_C_NONTARG", sgRNA) ~
          sub("CONTROL_C_([A-Za-z0-9]+)_.*", "\\1", sgRNA),
        grepl("CONTROL_C_NONTARG", sgRNA) ~ sgRNA,
        TRUE ~ NA_character_
      )
    )
  
  unique_entries <- unique(count_df_long_merged$unique_names)
  
  assignments <- sample(
    c(TRUE, FALSE),
    length(unique_entries),
    prob = c(percentage_kept_controls, percentage_used_controls),
    replace = TRUE
  )
  
  assignment_map <- data.frame(
    unique_names = unique_entries,
    keep = assignments
  )
  
  count_df_long_merged <- count_df_long_merged %>%
    dplyr::left_join(assignment_map, by = "unique_names")
  
  return_df <- count_df_long %>%
    dplyr::left_join(
      count_df_long_merged %>% dplyr::select(sgRNA, keep),
      by = "sgRNA"
    ) %>%
    dplyr::mutate(
      keep = ifelse(is.na(keep), FALSE, keep),
      group_category = ifelse(
        keep == TRUE,
        "kept_control",
        group_category
      )
    ) %>%
    dplyr::select(-keep)
  
  return(return_df)
}


subsample_controls_func <- function(count_df_long,
                                    cfg,
                                    merged_sgRNA_df = cfg$merged_sgRNA_df,
                                    percentage_kept_controls = 0.2) {
  
  percentage_used_controls <- 1 - percentage_kept_controls
  
  count_df_long_merged <- merge(
    count_df_long,
    merged_sgRNA_df[, c("sgrna_id", "entrez")],
    by.x = "sgRNA",
    by.y = "sgrna_id",
    all.x = TRUE
  )
  
  count_df_long_merged <- count_df_long_merged %>%
    dplyr::filter(group_category != "targeting")
  
  unique_entries <- unique(count_df_long_merged$sgRNA)
  
  assignments <- sample(
    c(TRUE, FALSE),
    length(unique_entries),
    prob = c(percentage_kept_controls, percentage_used_controls),
    replace = TRUE
  )
  
  assignment_map <- data.frame(
    sgRNA = unique_entries,
    keep = assignments
  )
  
  return_df <- count_df_long %>%
    dplyr::left_join(assignment_map, by = "sgRNA") %>%
    dplyr::mutate(
      keep = ifelse(is.na(keep), FALSE, keep),
      group_category = ifelse(
        keep == TRUE & group_category != "targeting",
        "kept_control",
        group_category
      )
    ) %>%
    dplyr::select(-keep)
  
  return(return_df)
}


count_df_long_to_wide <- function(count_df_long,
                                  print = TRUE,
                                  drop_0s = TRUE,
                                  recover_input = FALSE,
                                  for_sum = FALSE,
                                  ctrl_selection = c("targeting_control", "non_targeting_control")) {
  
  if (print == TRUE) {
    print(paste("NAs in count_df_long:", sum(is.na(count_df_long))))
  }
  
  maude_counts_df <- count_df_long %>%
    dplyr::mutate(
      isNontargeting = ifelse(group_category %in% ctrl_selection, TRUE, FALSE)
    ) %>%  
    dplyr::select(-c(group_category)) %>%
    tidyr::pivot_wider(names_from = condition, values_from = count) %>%
    dplyr::mutate(
      upper = ifelse(is.na(upper) & !is.na(lower), 0, upper),
      lower = ifelse(is.na(lower) & !is.na(upper), 0, lower)
    ) %>% 
    as.data.frame()
  
  if (print == TRUE) {
    print(paste("Dimension of wide_df:", paste(dim(maude_counts_df), collapse = " x ")))
    print(paste("NA's in maude_counts_df:", sum(is.na(maude_counts_df))))
    print(paste("NA's in Input:", sum(is.na(maude_counts_df$input))))
    print(paste("NA's in Upper:", sum(is.na(maude_counts_df$upper))))
    print(paste("NA's in Lower:", sum(is.na(maude_counts_df$lower))))
  } 
  
  if (recover_input == TRUE) {
    input_means <- maude_counts_df %>%
      dplyr::group_by(sgRNA, sublib) %>%
      dplyr::summarize(mean_input = mean(input, na.rm = TRUE), .groups = "drop")
    
    maude_counts_df <- maude_counts_df %>%
      dplyr::left_join(input_means, by = c("sgRNA", "sublib")) %>%
      dplyr::mutate(
        input = ifelse(
          is.na(input) & !is.na(upper) & !is.na(lower),
          mean_input,
          input
        )
      ) %>%
      dplyr::select(-mean_input)
    
  } else {
    if (for_sum == FALSE) {
      maude_counts_df <- maude_counts_df %>% 
        tidyr::drop_na()
    }
  }
  
  if (drop_0s == TRUE) {
    maude_counts_df <- maude_counts_df %>%
      dplyr::filter(!(input == 0 & upper == 0 & lower == 0))
  }
  
  if (for_sum == FALSE) {
    maude_counts_df <- maude_counts_df %>% 
      dplyr::mutate(
        input = input + 1,
        upper = upper + 1,
        lower = lower + 1
      )
  }
  
  if (print == TRUE) {
    print(paste("Dimension of wide_df:", paste(dim(maude_counts_df), collapse = " x ")))
    print(paste("NA's in maude_counts_df:", sum(is.na(maude_counts_df))))
    print(paste("NA's in Input:", sum(is.na(maude_counts_df$input))))
    print(paste("NA's in Upper:", sum(is.na(maude_counts_df$upper))))
    print(paste("NA's in Lower:", sum(is.na(maude_counts_df$lower))))
  }
  # MAUDE seems to not like floats, so we round to integers
  # TO DO: check how I can make maude use floats.
  
  maude_counts_df <- maude_counts_df %>% 
    dplyr::mutate(
      input = round(input),
      lower = round(lower),
      upper = round(upper)
    )
  
  return(maude_counts_df)
}