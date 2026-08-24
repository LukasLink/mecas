# R/R_functions/EA_32_make_waterfall_filename.R

make_waterfall_filename <- function(cfg) {
  
  parts <- c("waterfall")
  
  if (isTRUE(cfg$plots$waterfall$mark_cntrl)) {
    parts <- c(parts, "controls")
  }
  
  mark_special <- cfg$plots$waterfall$mark_special
  
  has_mark_special <- !is.null(mark_special) &&
    length(mark_special) > 0 &&
    !all(is.na(mark_special)) &&
    !all(mark_special == "")
  
  if (has_mark_special) {
    special <- paste(mark_special, collapse = "-")
    special <- gsub("[^A-Za-z0-9_.-]+", "_", special)
    parts <- c(parts, paste0("mark_", special))
  }
  
  if (cfg$plots$waterfall$mark_N_top_hits > 0) {
    parts <- c(parts, paste0("top", cfg$plots$waterfall$mark_N_top_hits))
  }
  
  if (!is.null(cfg$plots$waterfall$mark_all_signif_level) &&
      !is.na(cfg$plots$waterfall$mark_all_signif_level)) {
    parts <- c(parts, paste0("FDR", cfg$plots$waterfall$mark_all_signif_level))
  }
  
  if (isTRUE(cfg$plots$waterfall$no_text)) {
    parts <- c(parts, "no_text")
  }
  
  invisible(paste(parts, collapse = "_"))
}