# R/R_functions/EA_11_get_grouped_summary_wide.R

get_grouped_summary_wide <- function(count_df_long,
                                     stat = c("median", "mean")) {
  
  stat <- match.arg(stat)
  
  summary_fun <- switch(
    stat,
    median = function(x) median(x, na.rm = TRUE),
    mean = function(x) round(mean(x, na.rm = TRUE), 2)
  )
  
  count_df_long %>%
    dplyr::mutate(
      group_type = ifelse( group_category == "targeting", "targeting", "non-targeting")
    ) %>%
    dplyr::group_by(sublib, sample, group_type, bin_name) %>%
    dplyr::summarise(value = summary_fun(count), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = bin_name, values_from = value) %>%
    dplyr::mutate(group_id = paste(sublib, sample, group_type, sep = "_")) %>%
    dplyr::select(group_id, dplyr::everything(), -sublib, -sample, -group_type)
}
