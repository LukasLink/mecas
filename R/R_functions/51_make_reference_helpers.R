# 30_make_reference_helpers.R

write_star_fasta <- function(df, output_path) {
  
  required_cols <- c("sgrna_id", "full_oligo")
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop_log(
      "Cannot write STAR FASTA. Missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (any(is.na(df$sgrna_id) | df$sgrna_id == "")) {
    stop_log("Cannot write STAR FASTA. Some `sgrna_id` values are missing.")
  }
  
  if (any(is.na(df$full_oligo) | df$full_oligo == "")) {
    stop_log("Cannot write STAR FASTA. Some `full_oligo` values are missing.")
  }
  
  fasta_entries <- c(rbind(
    paste0(">", df$sgrna_id),
    df$full_oligo
  ))
  
  writeLines(fasta_entries, output_path)
}


build_star_gtf <- function(df) {
  
  required_cols <- c("sgrna_id", "full_oligo")
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop_log(
      "Cannot build STAR GTF. Missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  gtf_base <- data.frame(
    seqname = df$sgrna_id,
    source = "MECAS",
    start = 1L,
    end = nchar(df$full_oligo),
    score = ".",
    strand = "+",
    frame = ".",
    attribute = paste0(
      'gene_id "', df$sgrna_id, '"; ',
      'transcript_id "', df$sgrna_id, '"; ',
      'gene_name "', df$sgrna_id, '";'
    ),
    stringsAsFactors = FALSE
  )
  
  gene_df <- gtf_base
  gene_df$feature <- "gene"
  
  transcript_df <- gtf_base
  transcript_df$feature <- "transcript"
  
  exon_df <- gtf_base
  exon_df$feature <- "exon"
  
  gtf_df <- rbind(gene_df, transcript_df, exon_df)
  
  gtf_df <- gtf_df[
    order(gtf_df$seqname, match(gtf_df$feature, c("gene", "transcript", "exon"))),
  ]
  
  gtf_df <- gtf_df[
    ,
    c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  ]
  
  rownames(gtf_df) <- NULL
  
  gtf_df
}

slurm_mem_to_bytes <- function(mem) {
  mem <- trimws(as.character(mem))
  
  match <- regexec("^([0-9.]+)\\s*([KMGTP]?)(B)?$", mem, ignore.case = TRUE)
  parts <- regmatches(mem, match)[[1]]
  
  if (length(parts) == 0) {
    stop("Could not parse slurm_mem: ", mem)
  }
  
  value <- as.numeric(parts[2])
  unit <- toupper(parts[3])
  
  multiplier <- if (unit == "") {
    1
  } else {
    switch(
      unit,
      K = 1024,
      M = 1024^2,
      G = 1024^3,
      T = 1024^4,
      P = 1024^5,
      stop("Unsupported memory unit in slurm_mem: ", mem)
    )
  }
  
  floor(value * multiplier)
}