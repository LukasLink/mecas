# 31_run_make_reference.R

run_make_reference <- function(config_path, project_root_dir, cli_args) {
  #-----------------------------------------------------------------------------
  # Run setup
  #-----------------------------------------------------------------------------  
  config <- yaml::read_yaml(config_path)
  
  project_setup(
    project_root_dir = project_root_dir,
    config_path = config_path,
    setup_mode = "make_reference",
    use_old_suffix_construction = FALSE
  )
  
  use_pure_sgrna_sequence <- "--use-pure-sgrna-sequence" %in% cli_args
  stop_before_star <- "--only-ref-no-star" %in% cli_args
  
  #-----------------------------------------------------------------------------
  # Sanity Checks
  #----------------------------------------------------------------------------- 
  if (!exists("merged_sgRNA_df")) {
    stop_log("`merged_sgRNA_df` does not exist after project_setup(). Cannot make reference.")
  }
  
  required_cols <- c("sgrna_id", "seq", "align_seq")
  missing_cols <- setdiff(required_cols, colnames(merged_sgRNA_df))
  
  if (length(missing_cols) > 0) {
    stop_log(
      "merged_sgRNA_df is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!use_pure_sgrna_sequence && !"align_seq" %in% colnames(merged_sgRNA_df)) {
    stop_log(
      "`align_seq` column is missing from merged_sgRNA_df.\n",
      "Either provide `align_seq`, or rerun with `--use-pure-sgrna-sequence`."
    )
  }
  #-----------------------------------------------------------------------------
  # Use either align_seq or seq as the reference base
  #----------------------------------------------------------------------------- 
  ref_input_df <- merged_sgRNA_df
  
  if (use_pure_sgrna_sequence) {
    
    logger::log_warn(
      "Using pure sgRNA sequence from column `seq` for reference generation."
    )
    
    ref_input_df <- ref_input_df %>%
      mutate(full_oligo = seq)
    
  } else {
    
    n_missing_align_seq <- sum(is.na(ref_input_df$align_seq) | ref_input_df$align_seq == "")
    
    if (n_missing_align_seq > 0) {
      logger::log_warn(
        "Column `align_seq` contains {n_missing_align_seq} missing/empty values. ",
        "Using seq column to generate reference instead."
      )
      ref_input_df <- ref_input_df %>%
        mutate(full_oligo = align_seq)
      
    } else {
      ref_input_df <- ref_input_df %>%
        mutate(full_oligo = seq)      
    }
  }
  
  if (any(is.na(ref_input_df$full_oligo) | ref_input_df$full_oligo == "")) {
    stop_log(
      "Some rows still have missing/empty `full_oligo` after reference sequence selection."
    )
  }
  #-----------------------------------------------------------------------------
  # Generate Reference fasta and gtf files
  #----------------------------------------------------------------------------- 
  
  ref_gtf_path <- get_file_path(genome_output_folder, "ref.gtf")
  ref_fasta_path <- get_file_path(genome_output_folder, "ref.fa")
  
  
  write_star_fasta(ref_input_df, ref_fasta_path)
  
  star_gtf_df <- build_star_gtf(ref_input_df)
  
  write.table(
    star_gtf_df,
    ref_gtf_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  genome_expected_sjdb_overhang <- round(
    mean(nchar(ref_input_df$full_oligo)) - 1,
    0
  )
  
  total_ref_length <- sum(nchar(ref_input_df$full_oligo))
  n_refs <- nrow(ref_input_df)
  max_ref_length <- max(nchar(ref_input_df$full_oligo))
  
  genome_sa_index_n_bases <- round(
    min(14, log2(total_ref_length) / 2 - 1),
    0
  )
  
  genome_chr_bin_n_bits <- floor(
    min(18, log2(max(total_ref_length / n_refs, max_ref_length)))
  )
  
  genome_chr_bin_n_bits <- max(1, genome_chr_bin_n_bits)
  
  logger::log_info("Reference GTF written to: {ref_gtf_path}")
  logger::log_info("Reference FASTA written to: {ref_fasta_path}")
  logger::log_info("Suggested --genomeSAindexNbases: {genome_sa_index_n_bases}")
  logger::log_info("Suggested --sjdbOverhang: {genome_expected_sjdb_overhang}")
  logger::log_info("Suggested --genomeChrBinNbits: {genome_chr_bin_n_bits}")
  
  # Stop here if "--only-ref-no-star" flag is set
  if (stop_before_star){
    logger::log_info("Done.")
    return(invisible(list(
      ref_gtf_path = ref_gtf_path,
      ref_fasta_path = ref_fasta_path,
      genome_sa_index_n_bases = genome_sa_index_n_bases
    )))
  }
  logger::log_info("Proceeding to generate star_index with the suggested values.")
  #-----------------------------------------------------------------------------
  # Make slurm settings 
  #-----------------------------------------------------------------------------
  slurm_settings <- list(
    account    = slurm_account,
    qos        = slurm_qos,
    cpus       = slurm_cpus,
    mem        = slurm_mem,
    wall_time  = slurm_wall_time,
    partition  = slurm_partition,
    array      = 0,
    email      = slurm_email
  )
  #-----------------------------------------------------------------------------
  # Run STAR
  #-----------------------------------------------------------------------------
  star_index_output_dir <- make_clean_dir(opt$output_folder, "star_index")
  max_mem_for_star <- floor(slurm_mem_to_bytes(slurm_mem) * 0.90)
  
  run_shell_step(
    step_name = "generate_STAR_index",
    script_path = file.path(project_root_dir, "shell", "generate_STAR_index.sh"),
    args = c(
      "--runThreadN", slurm_cpus,                        # Goes directly to STAR
      "--runMode", "genomeGenerate",                     # Goes directly to STAR
      "--genomeDir", star_index_output_dir,              # Goes directly to STAR
      "--genomeFastaFiles", ref_fasta_path,              # Goes directly to STAR
      "--sjdbGTFfile", ref_gtf_path,                     # Goes directly to STAR
      "--limitGenomeGenerateRAM", max_mem_for_star,      # Goes directly to STAR
      "--genomeSAindexNbases", genome_sa_index_n_bases,  # Goes directly to STAR
      "--sjdbOverhang", genome_expected_sjdb_overhang,   # Goes directly to STAR
      "--genomeChrBinNbits", genome_chr_bin_n_bits,      # Goes directly to STAR
      "--use-modules", if (isTRUE(opt$use_modules)) "true" else "false",
      "--star-module", opt$star_module
    ),
    slurm_settings = slurm_settings,
    machine = opt$machine,
    log_dir = log_folder
  )
  logger::log_info("STAR index generated at: {star_index_output_dir}")
  logger::log_info("Copy this path to the config to proceed.")
  logger::log_info("Done.")
  return(invisible(list(
    ref_gtf_path = ref_gtf_path,
    ref_fasta_path = ref_fasta_path,
    star_index_output_dir = star_index_output_dir,
    genome_sa_index_n_bases = genome_sa_index_n_bases,
    genome_expected_sjdb_overhang = genome_expected_sjdb_overhang,
    max_mem_for_star = max_mem_for_star
  )))
}
