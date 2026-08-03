get_binStats <- function(bins, experiments, tail_p = 0.001) {
  
   if (!is.numeric(tail_p) || length(tail_p) != 1 || is.na(tail_p) || tail_p <= 0 || tail_p >= 0.5) {
    stop("`tail_p` must be one numeric value greater than 0 and smaller than 0.5.", call. = FALSE)
  }
  
  experiments <- unique(trimws(as.character(experiments)))
  experiments <- experiments[!is.na(experiments) & nzchar(experiments)]
  
  if (length(experiments) == 0) {
    stop("At least one experiment is required to construct `binStats`.", call. = FALSE)
  }
  
  binStats_base <- bins %>%
    dplyr::filter(sorted_or_unsorted == "sorted") %>%
    dplyr::arrange(bin_fraction_min, bin_fraction_max) %>%
    dplyr::mutate(
      fraction = bin_fraction_max - bin_fraction_min,
      binStartQ = dplyr::case_when(
        bin_fraction_min == 0 ~ tail_p,
        bin_fraction_max == 1 ~ 1 - tail_p - fraction,
        TRUE ~ bin_fraction_min
      ),
      binEndQ = dplyr::case_when(
        bin_fraction_min == 0 ~ tail_p + fraction,
        bin_fraction_max == 1 ~ 1 - tail_p,
        TRUE ~ bin_fraction_max
      )
    )
  
  invalid_adjusted_bounds <- binStats_base$binStartQ <= 0 |
    binStats_base$binEndQ >= 1 |
    binStats_base$binStartQ >= binStats_base$binEndQ
  
  if (any(invalid_adjusted_bounds)) {
    stop(
      "The MAUDE tail adjustment produced invalid quantile boundaries.\n",
      "Try using a smaller `tail_p` value.\n\n",
      "Affected bins:\n",
      paste0("  - ", binStats_base$bin_name[invalid_adjusted_bounds], collapse = "\n"),
      call. = FALSE
    )
  }
  
  binStats_base <- binStats_base %>%
    dplyr::transmute(
      Bin = bin_name,
      binStartQ = binStartQ,
      binEndQ = binEndQ,
      fraction = fraction,
      binStartZ = stats::qnorm(binStartQ),
      binEndZ = stats::qnorm(binEndQ)
    ) %>%
    as.data.frame()
  
  binStats <- dplyr::bind_rows(
    lapply(
      experiments,
      function(current_exp) {
        binStats_base %>%
          dplyr::mutate(exp = current_exp)
      }
    )
  )
  
  binStats <- binStats %>%
    dplyr::select(Bin, binStartQ, binEndQ, fraction, binStartZ, binEndZ, exp) %>%
    as.data.frame()
  
  rownames(binStats) <- NULL
  
  binStats
}
