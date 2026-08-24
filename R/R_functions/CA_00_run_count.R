# R/R_functions/CA_00_run_count.R

run_count <- function(cfg) {
  
  
  if(isTRUE(cfg$counting$umis_as_sublibs)){
    log_info("Starting Consolidation of UMI counts...")
    count_df_long <- get_reads_per_umi_count_df_long_and_make_coverage_file(cfg = cfg)
  } else {
    log_info("Starting conversion of Read/UMI counts from bcwithqc sparse matrices...")
    count_df_long <- get_count_df_long_and_make_coverage_file(cfg = cfg)
  }
  
  # If we are processing umis, might as well already provide the reads_counts too. 
  if (identical("umis", cfg$counting$data_type)){
    count_df_long_reads <- get_read_count_df_long(cfg, return_df = FALSE)
  }
  
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