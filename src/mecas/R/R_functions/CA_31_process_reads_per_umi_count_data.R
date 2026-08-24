# R/R_functions/CA_31_process_reads_per_umi_count_data.R

process_reads_per_umi_count_data <- function(
    cfg,
    manifest = cfg$manifest
) {
  
  pipeline_names <- unique(manifest$pipeline_name)
  count_df_long <- NULL
  
  for (pipeline_name in pipeline_names) {
    
    rds_path <- file.path(cfg$paths$reads_per_umi_count, paste0(pipeline_name, ".rds"))
    
    if (!file.exists(rds_path)) {
      stop_log(
        "Expected reads-per-UMI count file does not exist for `",
        pipeline_name,
        "`:\n  ",
        rds_path,
        call. = FALSE
      )
    }
    
    df <- tryCatch(
      readRDS(rds_path),
      error = function(e) {
        stop_log(
          "Failed to import reads-per-UMI counts for `",
          pipeline_name,
          "`:\n",
          conditionMessage(e)
        )
      }
    )
    # --------------------------------------------------------------------------
    # Safety Checks
    # --------------------------------------------------------------------------
    if (!is.data.frame(df)) {
      stop_log("The following RDS does not contain a data frame:\n", rds_path)
    }
    
    required_columns <- c("sgRNA", "sublib", "bin_name", "count", "group_category", "sample", "exp")
    missing_columns <- setdiff(required_columns, colnames(df))
    
    if (length(missing_columns) > 0) {
      stop_log(
        "The RDS file:\n", rds_path, "\n", "is missing required columns:\n",
        paste(missing_columns, collapse = ", ")
      )
    }
    if (nrow(df) == 0){
      log_warn("The file:\n{rds_path}\n contains no data.")
    }
    if (nrow(df) < 100){
      log_warn("The file:\n{rds_path}\n contains less than 100 entries.")
    }
    # --------------------------------------------------------------------------   
    if(is.null(count_df_long)){
      # First itteration create count_df_long
      count_df_long <- df
    } else {
      # Every other itteration: add rows to count_df_long
      count_df_long <- dplyr::bind_rows(count_df_long, df)
    }
  } # end of iterating through pipeline_names
  
  return(count_df_long)
}