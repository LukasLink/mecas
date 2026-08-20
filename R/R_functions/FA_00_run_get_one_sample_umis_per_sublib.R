# R/R_functions/FA_00_run_get_one_sample_umis_per_sublib.R


run_get_one_sample_umis_per_sublib <- function(cfg, sample){
  
  manifest <- load_existing_tsv_or_stop(
    cfg$paths$manifest,
    "standardized FASTQ manifest"
  )
  
  if(!sample %in% manifest$pipeline_name){
    stop_log("In run_get_one_sample_umis_per_sublib: ",sample, " was not found in the manifest")
  }
  
  manifest_row <- manifest[manifest$pipeline_name == sample, ][1, ]
  
  bin_name <- manifest_row$bin_name
  sample_name <- manifest_row$sample
  
  bam_fpath <- file.path(
    cfg$paths$bcwithqc_output_folder,
    sample,
    "with_bc_umi.sorted.bam"
  )
  
  reads_per_umi_count_df <- build_reads_per_umi_count_df(
    bam_fpath = bam_fpath,
    library = cfg$merged_sgRNA_df,
    set_bin_name_to = bin_name,
    set_sample_to = sample_name
  )
  
  dir.create(cfg$paths$snake$count_df_meta, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$paths$reads_per_umi_count, recursive = TRUE, showWarnings = FALSE)
  
  writeLines(
    as.character(nrow(reads_per_umi_count_df)),
    file.path(cfg$paths$snake$count_df_meta, paste0(sample, ".txt"))
  )
  
  saveRDS(reads_per_umi_count_df, file.path(cfg$paths$reads_per_umi_count, paste0(sample,".rds")))
}