# This function is deprecated and it should be tested if it can just be removed without breaking anything. 

make_star_custom_reference <- function(merged_sgRNA_df,
                                       genome_output_folder) {
  
  fasta_path <- file.path(genome_output_folder, "reference.fa")
  gtf_path <- file.path(genome_output_folder, "reference.gtf")
  
  # ---------------------------------------------------------------------------
  # Sanity checks
  # ---------------------------------------------------------------------------
  
  required_cols <- c("sgrna_id", "align_seq")
  missing_cols <- setdiff(required_cols, colnames(merged_sgRNA_df))
  
  if (length(missing_cols) > 0) {
    stop_log(
      "Cannot make custom STAR reference. Missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (any(is.na(merged_sgRNA_df$sgrna_id) | merged_sgRNA_df$sgrna_id == "")) {
    stop_log("Cannot make custom STAR reference. Some `sgrna_id` values are missing.")
  }
  
  if (any(is.na(merged_sgRNA_df$align_seq) | merged_sgRNA_df$align_seq == "")) {
    stop_log(
      "A column with alignment sequences must be included in the library file for make_reference. ",
      "(see library_formatting_requirements.txt)"
    )
  }
  
  # ---------------------------------------------------------------------------
  # Prepare reference data
  # ---------------------------------------------------------------------------
  
  ref_df <- merged_sgRNA_df %>%
    dplyr::transmute(
      sgrna_id = sgrna_id,
      full_oligo = toupper(align_seq)
    )
  
  # ---------------------------------------------------------------------------
  # Write FASTA/GTF using shared helpers
  # ---------------------------------------------------------------------------
  
  write_star_fasta(ref_df, fasta_path)
  
  gtf_df <- build_star_gtf(ref_df)
  
  write.table(
    gtf_df,
    gtf_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # Calculate STAR genomeSAindexNbases
  # ---------------------------------------------------------------------------
  
  total_reference_length <- sum(nchar(ref_df$full_oligo))
  
  genome_sa_index_n_bases <- round(
    min(14, log2(total_reference_length) / 2 - 1),
    0
  )
  
  genome_sa_index_n_bases <- max(1, as.integer(genome_sa_index_n_bases))
  
  logger::log_info("Wrote FASTA: {fasta_path}")
  logger::log_info("Wrote GTF: {gtf_path}")
  logger::log_info("For STAR --genomeSAindexNbases use: {genome_sa_index_n_bases}")
  
  invisible(list(
    fasta_path = fasta_path,
    gtf_path = gtf_path,
    genome_sa_index_n_bases = genome_sa_index_n_bases
  ))
}