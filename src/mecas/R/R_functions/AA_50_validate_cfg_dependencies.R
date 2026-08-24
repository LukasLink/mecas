# R/R_functions/AA_50_validate_cfg_dependencies.R

validate_cfg_dependencies <- function(cfg){
  if (isTRUE(cfg$counting$umis_as_sublibs)){
    if (identical("reads", cfg$counting$data_type)){
      stop_log(
        "data_type = 'reads' is not supported in '--umis-as-sublibs' mode.",
        "Please use '--reads-or-umis umis' instead.")
    }
    if (identical("rep_sublib", cfg$replicates$method)){
      stop_log(
        "Replicate handeling method 'rep_sublib' is not supported in '--umis-as-sublibs' mode.")
    }
    
    if (
      identical("sublib", cfg$replicates$combine_for_guide_stats)
      |
      identical("sublib", cfg$replicates$combine_for_gene_stats)
      ){
      log_warn(
        "WARNING! Either combine_for_guide_stats or combine_for_gene_stats was set to 'sublib' in '--umis-as-sublibs' mode. ",
        "This might lead to unintended behaviour or nullify the point of using that mode to begin with.",
        "It SHOULD not break anything, but don't use this combination unless you have good reason!")
    }
  }
}