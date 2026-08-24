# R/R_functions/DA_12_subsample_controls_func.R

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