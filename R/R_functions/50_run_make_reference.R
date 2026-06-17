# 31_run_make_reference.R

run_make_reference <- function(config_path, project_root_dir, cli_args, run_setup = TRUE, cfg = NULL) {
  #-----------------------------------------------------------------------------
  # Run setup
  #-----------------------------------------------------------------------------
  if (isTRUE(run_setup)){
    cfg <- project_setup(
      project_root_dir = project_root_dir,
      config_path = config_path,
      setup_mode = "make_reference",
      use_old_suffix_construction = FALSE
    )
  } else {
    cfg$slurm$array <- 1
  }
  if (is.null(cfg)){
    stop("Config is NULL after run_make_reference stup. This is a code Error, users can not fix this.")
  }
  
  # This is the old version from the external function, the internal should always use just the sgRNA sequence. 
  # use_pure_sgrna_sequence <- "--use-pure-sgrna-sequence" %in% cli_args
  use_pure_sgrna_sequence <- TRUE
  stop_before_star <- "--only-ref-no-star" %in% cli_args
  
  #-----------------------------------------------------------------------------
  # Sanity Checks
  #-----------------------------------------------------------------------------
  
  if (is.null(cfg$merged_sgRNA_df)) {
    stop_log("`cfg$merged_sgRNA_df` does not exist after project_setup(). Cannot make reference.")
  }
  
  required_cols <- c("sgrna_id", "seq", "align_seq")
  missing_cols <- setdiff(required_cols, colnames(cfg$merged_sgRNA_df))
  
  if (length(missing_cols) > 0) {
    stop_log(
      "cfg$merged_sgRNA_df is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!use_pure_sgrna_sequence && !"align_seq" %in% colnames(cfg$merged_sgRNA_df)) {
    stop_log(
      "`align_seq` column is missing from cfg$merged_sgRNA_df.\n",
      "Either provide `align_seq`, or rerun with `--use-pure-sgrna-sequence`."
    )
  }
  
  #-----------------------------------------------------------------------------
  # Use either align_seq or seq as the reference base
  #-----------------------------------------------------------------------------
  
  ref_input_df <- cfg$merged_sgRNA_df
  
  if (use_pure_sgrna_sequence) {
    
    logger::log_warn(
      "Using pure sgRNA sequence from column `seq` for reference generation."
    )
    
    ref_input_df <- ref_input_df %>%
      dplyr::mutate(full_oligo = seq)
    
  } else {
    
    n_missing_align_seq <- sum(
      is.na(ref_input_df$align_seq) | ref_input_df$align_seq == ""
    )
    
    if (n_missing_align_seq > 0) {
      logger::log_warn(
        "Column `align_seq` contains {n_missing_align_seq} missing/empty values. ",
        "Using `seq` column to generate reference instead."
      )
      
      ref_input_df <- ref_input_df %>%
        dplyr::mutate(full_oligo = seq)
      
    } else {
      ref_input_df <- ref_input_df %>%
        dplyr::mutate(full_oligo = align_seq)
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
  
  ref_gtf_path <- get_file_path(cfg$paths$genome_output_folder, "ref.gtf")
  ref_fasta_path <- get_file_path(cfg$paths$genome_output_folder, "ref.fa")
  
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
  
  genome_sa_index_n_bases <- floor(
    min(14, log2(total_ref_length) / 2 - 1)
  )
  
  genome_sa_index_n_bases <- max(1, genome_sa_index_n_bases)
  
  genome_chr_bin_n_bits <- floor(
    min(18, log2(max(total_ref_length / n_refs, max_ref_length)))
  )
  
  genome_chr_bin_n_bits <- max(1, genome_chr_bin_n_bits)
  
  logger::log_info("Reference GTF written to: {ref_gtf_path}")
  logger::log_info("Reference FASTA written to: {ref_fasta_path}")
  logger::log_info("Suggested --genomeSAindexNbases: {genome_sa_index_n_bases}")
  logger::log_info("Suggested --sjdbOverhang: {genome_expected_sjdb_overhang}")
  logger::log_info("Suggested --genomeChrBinNbits: {genome_chr_bin_n_bits}")
  
  #-----------------------------------------------------------------------------
  # Stop here if "--only-ref-no-star" flag is set
  #-----------------------------------------------------------------------------
  
  if (isTRUE(stop_before_star)) {
    logger::log_info("Done.")
    
    return(invisible(list(
      cfg = cfg,
      ref_gtf_path = ref_gtf_path,
      ref_fasta_path = ref_fasta_path,
      genome_sa_index_n_bases = genome_sa_index_n_bases,
      genome_expected_sjdb_overhang = genome_expected_sjdb_overhang,
      genome_chr_bin_n_bits = genome_chr_bin_n_bits
    )))
  }
  
  logger::log_info("Proceeding to generate STAR index with the suggested values.")
  
  #-----------------------------------------------------------------------------
  # Run STAR
  #-----------------------------------------------------------------------------
  
  star_index_output_dir <- make_clean_dir(
    cfg$paths$output_folder,
    "star_index"
  )
  
  max_mem_for_star <- floor(
    slurm_mem_to_bytes(cfg$slurm$mem) * 0.90
  )
  
  run_shell_step(
    step_name = "generate_STAR_index",
    script_path = file.path(project_root_dir, "shell", "generate_STAR_index.sh"),
    cfg = cfg,
    args = c(
      "--runThreadN", as.character(cfg$slurm$cpus),
      "--runMode", "genomeGenerate",
      "--genomeDir", star_index_output_dir,
      "--genomeFastaFiles", ref_fasta_path,
      "--sjdbGTFfile", ref_gtf_path,
      "--limitGenomeGenerateRAM", as.character(max_mem_for_star),
      "--genomeSAindexNbases", as.character(genome_sa_index_n_bases),
      "--sjdbOverhang", as.character(genome_expected_sjdb_overhang),
      "--genomeChrBinNbits", as.character(genome_chr_bin_n_bits),
      "--use-modules", if (isTRUE(cfg$modules$use_modules)) "true" else "false",
      "--star-module", cfg$modules$star
    ),
    log_dir = cfg$paths$log_folder %||% cfg$paths$logs_folder
  )
  
  logger::log_info("STAR index generated at: {star_index_output_dir}")
  if (isTRUE(run_setup)){
    logger::log_info("Copy this path to the config to proceed.")
    logger::log_info("Done.")    
  }
  
  return(invisible(list(
    cfg = cfg,
    ref_gtf_path = ref_gtf_path,
    ref_fasta_path = ref_fasta_path,
    star_index_output_dir = star_index_output_dir,
    genome_sa_index_n_bases = genome_sa_index_n_bases,
    genome_expected_sjdb_overhang = genome_expected_sjdb_overhang,
    genome_chr_bin_n_bits = genome_chr_bin_n_bits,
    max_mem_for_star = max_mem_for_star
  )))
}