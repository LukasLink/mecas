# R/R_functions/CA_00_run_count.R

run_count <- function(cfg) {
  
  log_info("Starting conversion of Read/UMI counts from bcwithqc sparse matrix...")
  
  count_df_long <- get_count_df_long_and_make_coverage_file(cfg = cfg)
  
  count_df_tsv_gz <- file.path(
    cfg$paths$results_output_folder,
    "count_df.tsv.gz"
  )
  
  # R-specific compact version
  saveRDS(
    count_df_long,
    file = cfg$paths$count_df_fpath,
    compress = "gzip"
  )
  
  # Compressed, portable TSV version
  readr::write_tsv(
    count_df_long,
    file = count_df_tsv_gz
  )
  
  log_info(
    "Read/UMI counts written to {cfg$paths$count_df_fpath} and {count_df_tsv_gz}"
  )
}