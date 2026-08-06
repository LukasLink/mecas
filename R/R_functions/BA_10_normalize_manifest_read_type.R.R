# R/R_functions/BA_10_normalize_manifest_read_type.R


normalize_manifest_read_type <- function(manifest, read_col = "read") {
  # Backwards compatibility with older single-end manifests
  if (!read_col %in% colnames(manifest)) {
    return(rep("SE", nrow(manifest)))
  }
  
  read_type <- toupper(trimws(as.character(manifest[[read_col]])))
  
  # Empty or NA means single-end
  read_type[is.na(read_type) | read_type == ""] <- "SE"
  
  invalid_read_types <- setdiff(unique(read_type), c("SE", "R1", "R2"))
  
  if (length(invalid_read_types) > 0) {
    stop_log(
      "Manifest contains invalid values in column `",
      read_col,
      "`: ",
      paste(invalid_read_types, collapse = ", "),
      ". Expected empty values, R1, or R2."
    )
  }
  
  read_type
}


