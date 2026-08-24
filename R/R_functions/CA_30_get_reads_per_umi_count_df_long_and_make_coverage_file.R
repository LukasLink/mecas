# R/R_functions/CA_30_get_reads_per_umi_count_df_long_and_make_coverage_file.R

get_reads_per_umi_count_df_long_and_make_coverage_file <- function(cfg) {
  logger::log_info("Reading reads-per-UMI counts...")
  
  count_df_long <- process_reads_per_umi_count_data(cfg = cfg)
  
  logger::log_info("Finished reading reads-per-UMI counts.")
  
  logger::log_info("Creating reads-per-UMI coverage file...")
  
  coverage_df <- make_reads_per_umi_coverage_df(count_df_long, cfg = cfg)
  
  writexl::write_xlsx(
    coverage_df,
    file.path(
      cfg$paths$results_output_folder,
      "mapping_results.xlsx"
    )
  )
  
  logger::log_info("Finished creating reads-per-UMI coverage file.")
  logger::log_info(
    paste(
      "The coverage file is:",
      file.path(
        cfg$paths$results_output_folder,
        "mapping_results.xlsx"
      )
    )
  )
  
  count_df_long <- apply_optional_count_df_operations(count_df_long, cfg = cfg)
  
  return(count_df_long)
}