# R/R_functions/DA_11_check_all_bins_present.R

check_all_bins_present <- function(count_df_long, cfg) {
  
  required_count_columns <- c("sublib", "sample", "bin_name")
  missing_count_columns <- base::setdiff(required_count_columns, colnames(count_df_long))
  
  if (length(missing_count_columns) > 0) {
    stop_log(
      "`count_df_long` is missing columns required to check MAUDE bins:\n",
      paste0("  - ", missing_count_columns, collapse = "\n")
    )
  }
  
  if (!is.data.frame(cfg$bins) || !"bin_name" %in% colnames(cfg$bins)) {
    stop_log(
      "`cfg$bins` must be a validated data frame containing a `bin_name` column."
    )
  }
  
  required_bins <- unique(as.character(cfg$bins$bin_name))
  
  if (length(required_bins) == 0 || any(is.na(required_bins) | !nzchar(required_bins))) {
    stop_log("No valid bin names were found in `cfg$bins`.")
  }
  
  observed_bin_names <- unique(stats::na.omit(as.character(count_df_long$bin_name)))
  unexpected_bins <- base::setdiff(observed_bin_names, required_bins)
  
  if (length(unexpected_bins) > 0) {
    stop_log(
      "`count_df_long` contains bin names that are not defined in `cfg$bins`:\n",
      paste0("  - ", unexpected_bins, collapse = "\n")
    )
  }
  
  missing_bins_df <- count_df_long %>%
    dplyr::distinct(sublib, sample, bin_name) %>%
    dplyr::group_by(sublib, sample) %>%
    dplyr::summarise(
      present_bins = list(unique(stats::na.omit(as.character(bin_name)))),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      missing_bins = purrr::map(
        present_bins,
        ~ base::setdiff(required_bins, .x)
      )
    ) %>%
    dplyr::filter(lengths(missing_bins) > 0)
  
  if (nrow(missing_bins_df) > 0) {
    missing_details <- missing_bins_df %>%
      dplyr::mutate(
        missing_bins_text = purrr::map_chr(
          missing_bins,
          ~ paste(.x, collapse = ", ")
        ),
        message = paste0(
          "  - sample='", sample,
          "', sublib='", sublib,
          "': missing ", missing_bins_text
        )
      ) %>%
      dplyr::pull(message)
    
    stop_log(
      "MAUDE requires all defined bins to be provided for every sample and sublibrary.\n",
      "Required bins: ",
      paste(required_bins, collapse = ", "),
      "\n\nMissing bins by sample and sublibrary:\n",
      paste(missing_details, collapse = "\n")
    )
  }
  
  invisible(TRUE)
}