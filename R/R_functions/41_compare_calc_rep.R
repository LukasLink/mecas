# 45_compare_calc_rep.R
#-----------------------------------------------------------------------------
# find_Hits_and_FDR_appearing_in_X_replicates
#-----------------------------------------------------------------------------
find_Hits_and_FDR_appearing_in_X_replicates <- function(df_list,
                                                        name_list,
                                                        FDR_threshold = 0.05,
                                                        X = NULL){
  if (is.null(X)){
    X <- length(df_list)
  }
  filtered_list <- vector("list", length(df_list))
  for (i in seq_along(df_list)){
    filtered_list[[i]] <- df_list[[i]] %>%
      apply_FDR_filter(FDR_threshold) %>% 
      mutate(source = name_list[i])
  }
  
  # Combine
  combined <- dplyr::bind_rows(filtered_list)
  
  # Count in how many *sources* each entrez appears (after filtering)
  entrez_keep <- combined %>%
    dplyr::distinct(entrez, source) %>%
    dplyr::count(entrez, name = "n_sources") %>%
    dplyr::filter(n_sources >= X) %>%
    dplyr::pull(entrez)
  
  # For kept entrez: choose the row with the highest FDR (ties broken deterministically)
  results_df <- combined %>%
    dplyr::filter(entrez %in% entrez_keep) %>%
    dplyr::arrange(entrez, dplyr::desc(FDR), source) %>%
    dplyr::group_by(entrez) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  new_FDR <- if (nrow(results_df) == 0) NA_real_ else max(results_df$FDR, na.rm = TRUE)
  return(list(results_df = results_df, new_FDR = new_FDR))
}

automate_calc_replicate_comparison <- function(cfg,
                                               file_suffix = cfg$suffix$file_suffix,
                                               rds_output_folder = cfg$paths$rds_output_folder,
                                               FDR_threshold = 0.05,
                                               hits_in_X_reps = NULL,
                                               correlation_heatmap = TRUE,
                                               overlap = TRUE,
                                               venn_diagram = TRUE) {
  
  Hits_rep_0 <- add_info_wrapper(
    file_suffix = file_suffix,
    cfg = cfg,
    folder = rds_output_folder
  )
  
  # Not using the normal file_info_suffix from cfg here so it can be overwritten
  file_info_suffix <- get_file_info_suffix(file_suffix)
  
  # List all files in the rds output folder
  rds_files <- list.files(
    rds_output_folder,
    pattern = paste0(".*", file_info_suffix, "_rep[0-9]+\\.rds$"),
    full.names = TRUE
  )
  
  if (length(rds_files) == 0) {
    stop_log(
      "No consensus replicate RDS files found in:\n  ",
      rds_output_folder,
      "\n\nExpected files matching:\n  ",
      paste0(".*", file_info_suffix, "_rep[0-9]+\\.rds$")
    )
  }
  
  # Extract everything between file_info_suffix + "_" and ".rds"
  # e.g. "...<file_info_suffix>_no_sublib_L1.rds" -> "no_sublib_L1"
  rep_name <- sub(
    paste0(".*", file_info_suffix, "_rep(.+)\\.rds$"),
    "\\1",
    basename(rds_files)
  )
  
  rep_name <- unique(rep_name)
  rep_name_long <- paste0("calc_rep_", rep_name)
  
  if (file_info_suffix == ""){
    rep_file_suffixes <- sub(
      paste0("^.*(", file_info_suffix, "rep.+\\.rds$)"),
      "_\\1",
      basename(rds_files)
    )
  } else {
    rep_file_suffixes <- sub(
      paste0("^.*(", file_info_suffix, "_rep.+\\.rds$)"),
      "_\\1",
      basename(rds_files)
    )
  }

  
  # Generate Hits_df_list using add_info_wrapper for each sublib_suffixes
  Hits_df_list <- lapply(rep_file_suffixes, function(rep_file_suffix) {
    add_info_wrapper(
      file_suffix = rep_file_suffix,
      cfg = cfg,
      folder = rds_output_folder
    )
  })
  
  Hits_and_FDR <- find_Hits_and_FDR_appearing_in_X_replicates(
    df_list = c(list(Hits_rep_0), Hits_df_list),
    name_list = c("calc_rep_0", rep_name_long),
    X = hits_in_X_reps,
    FDR_threshold = FDR_threshold
  )
  
  return(list(
    Hits_in_X_df = Hits_and_FDR$results_df,
    new_FDR = Hits_and_FDR$new_FDR
  ))
}

#-----------------------------------------------------------------------------
# automate_calc_replicate_means
#-----------------------------------------------------------------------------
automate_calc_replicate_means <- function(cfg,
                                          file_suffix = cfg$suffix$file_suffix,
                                          rds_output_folder = cfg$paths$rds_output_folder) {
  
  Hits_rep_0 <- add_info_wrapper(
    file_suffix = file_suffix,
    cfg = cfg,
    folder = rds_output_folder
  )
  
  # Not using the normal file_info_suffix from cfg here so it can be overwritten
  file_info_suffix <- get_file_info_suffix(file_suffix)
  
  rds_files <- list.files(
    rds_output_folder,
    pattern = paste0(".*", file_info_suffix, "_rep[0-9]+\\.rds$"),
    full.names = TRUE
  )
  
  if (length(rds_files) == 0) {
    stop_log(
      "No consensus replicate RDS files found in:\n  ",
      rds_output_folder,
      "\n\nExpected files matching:\n  ",
      paste0(".*", file_info_suffix, "_rep[0-9]+\\.rds$")
    )
  }
  
  rep_name <- sub(
    paste0(".*", file_info_suffix, "_rep(.+)\\.rds$"),
    "\\1",
    basename(rds_files)
  )
  
  rep_name <- unique(rep_name)
  
  if (file_info_suffix == ""){
    rep_file_suffixes <- sub(
      paste0("^.*(", file_info_suffix, "rep.+\\.rds$)"),
      "_\\1",
      basename(rds_files)
    )
  } else {
    rep_file_suffixes <- sub(
      paste0("^.*(", file_info_suffix, "_rep.+\\.rds$)"),
      "_\\1",
      basename(rds_files)
    )
  }
  
  Hits_df_list <- lapply(rep_file_suffixes, function(rep_file_suffix) {
    logger::log_info("Processing replicate suffix: {rep_file_suffix}")
    
    add_info_wrapper(
      file_suffix = rep_file_suffix,
      cfg = cfg,
      folder = rds_output_folder
    )
  })
  
  all_hits_df <- dplyr::bind_rows(c(list(Hits_rep_0), Hits_df_list))
  
  mean_cols <- c("stoufferZ", "meanZ", "significanceZ", "p.value", "FDR")
  mean_cols <- base::intersect(mean_cols, colnames(all_hits_df))
  
  other_cols <- base::setdiff(colnames(all_hits_df), c("entrez", mean_cols))
  
  Mean_Hits_df <- all_hits_df %>%
    dplyr::group_by(entrez) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(mean_cols),
        ~ mean(.x, na.rm = TRUE)
      ),
      dplyr::across(
        dplyr::all_of(other_cols),
        ~ dplyr::first(.x)
      ),
      .groups = "drop"
    ) %>%
    dplyr::select(dplyr::all_of(c("entrez", other_cols, mean_cols)))
  
  return(Mean_Hits_df)
}