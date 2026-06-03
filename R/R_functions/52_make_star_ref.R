make_star_custom_reference <- function(merged_sgRNA_df = get("merged_sgRNA_df", envir = .GlobalEnv)) {

  genome_output_folder <- get("genome_output_folder", envir = .GlobalEnv)
  
  fasta_path <- file.path(genome_output_folder, "reference.fa")
  gtf_path   <- file.path(genome_output_folder, "reference.gtf")
  
  # ---------------------------------------------------------------------------
  # Sanity checks
  # ---------------------------------------------------------------------------

  if (any(is.na(merged_sgRNA_df$align_seq)) || any(merged_sgRNA_df$align_seq == "")) {
    stop("A column with alignment sequences must be included in the library file for make_reference. (see library_formatting_requirements.txt)")
  }
  
  # Keep only the required reference information
  ref_df <- merged_sgRNA_df %>%
    dplyr::mutate(
      align_seq = toupper(align_seq),
      seq_length = nchar(align_seq)
    )
  
  # ---------------------------------------------------------------------------
  # Write FASTA
  # ---------------------------------------------------------------------------
  fasta_entries <- c(rbind(
    paste0(">", ref_df$sgrna_id),
    ref_df$align_seq
  ))
  
  writeLines(fasta_entries, fasta_path)
  
  # ---------------------------------------------------------------------------
  # Write minimal GTF
  #
  # STAR only needs valid GTF structure and the core attributes gene_id and
  # transcript_id. Here each sgRNA/reference sequence is treated as one exon.
  # Using sgrna_id for gene_id prevents multiple sgRNAs targeting the same gene
  # from being collapsed at the annotation level.
  # ---------------------------------------------------------------------------
  gtf_df <- ref_df %>%
    dplyr::transmute(
      seqname = sgrna_id,
      source = "custom_reference",
      feature = "exon",
      start = 1L,
      end = seq_length,
      score = ".",
      strand = "+",
      frame = ".",
      attribute = paste0(
        'gene_id "', sgrna_id, '"; ',
        'transcript_id "', sgrna_id, '";'
      )
    )
  
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
  total_reference_length <- sum(ref_df$seq_length)
  
  genomeSAindexNbases <- round(
    min(14, log2(total_reference_length) / 2 - 1),
    0
  )
  
  genomeSAindexNbases <- max(1, as.integer(genomeSAindexNbases))
  
  message("Wrote FASTA: ", fasta_path)
  message("Wrote GTF:   ", gtf_path)
  message("For STAR --genomeSAindexNbases use: ", genomeSAindexNbases)
  
  return(genomeSAindexNbases)
}