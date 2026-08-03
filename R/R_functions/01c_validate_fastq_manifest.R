# One manifest row represents one FASTQ file.
# Single-end sample: one row with read = NA.
# Paired-end sample: two rows with the same pipeline_name and read = R1/R2.

validate_fastq_manifest <- function(manifest) {
  
  required_cols <- c(
    "bin_name",
    "sublibrary",
    "sample",
    "pipeline_name",
    "read",
    "fastq_id"
  )
  
  missing_cols <- setdiff(required_cols, colnames(manifest))
  
  if (length(missing_cols) > 0) {
    stop_log(
      "Manifest is missing columns required to validate FASTQ grouping:\n",
      paste0("  - ", missing_cols, collapse = "\n")
    )
  }
  
  groups <- split(seq_len(nrow(manifest)), manifest$pipeline_name)
  invalid_groups <- character(0)
  
  for (pipeline_name in names(groups)) {
    idx <- groups[[pipeline_name]]
    reads <- manifest$read[idx]
    
    valid_single_end <-
      length(idx) == 1 &&
      (is.na(reads[1]) || !nzchar(reads[1]))
    
    valid_paired_end <-
      length(idx) == 2 &&
      !any(is.na(reads)) &&
      !anyDuplicated(reads) &&
      setequal(reads, c("R1", "R2"))
    
    if (!valid_single_end && !valid_paired_end) {
      read_description <- paste(
        ifelse(is.na(reads) | !nzchar(reads), "<single-end>", reads),
        collapse = ", "
      )
      
      invalid_groups <- c(
        invalid_groups,
        paste0(
          pipeline_name,
          ": ",
          length(idx),
          " file(s); read values = ",
          read_description
        )
      )
    }
  }
  
  if (length(invalid_groups) > 0) {
    stop_log(
      "Invalid FASTQ grouping in the manifest. Each pipeline sample must have either:\n",
      "  - exactly one row with an empty `read` value for single-end data, or\n",
      "  - exactly two rows with `read` values R1 and R2 for paired-end data.\n\n",
      paste0("  - ", invalid_groups, collapse = "\n")
    )
  }
  
  if (anyDuplicated(manifest$fastq_id)) {
    duplicated_ids <- unique(manifest$fastq_id[duplicated(manifest$fastq_id)])
    
    stop_log(
      "Manifest contains duplicate `fastq_id` values:\n",
      paste0("  - ", duplicated_ids, collapse = "\n")
    )
  }
  
  invisible(TRUE)
}