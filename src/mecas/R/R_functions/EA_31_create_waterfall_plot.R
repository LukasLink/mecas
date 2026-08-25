# R/R_functions/EA_31_create_waterfall_plot.R

create_waterfall_plot <- function(Hits_df, cfg) {
  
  include_controls_list <- cfg$controls$include_controls %||% character()
  
  if (length(include_controls_list) > 0) {
    for (control_gene in include_controls_list) {
      control_gene <- sub("_[^_]*$", "", control_gene)
      Hits_df$symbol[
        grepl(control_gene, Hits_df$entrez)
      ] <- control_gene
    }
  }
  
  p <- plot_significance_by_rank_from_cfg(Hits_df, cfg)

  invisible(p)
}

plot_significance_by_rank_from_cfg <- function(Hits_df, cfg) {
  p <- plot_significance_by_rank(
    Hits_df = Hits_df,
    mark_cntrl = cfg$plots$waterfall$mark_cntrl,
    mark_special = cfg$plots$waterfall$mark_special,
    mark_N_top_hits = cfg$plots$waterfall$mark_N_top_hits,
    box_padding = cfg$plots$waterfall$box_padding,
    no_text = cfg$plots$waterfall$no_text,
    signif_lines = cfg$plots$waterfall$signif_lines,
    mark_all_signif_level = cfg$plots$waterfall$mark_all_signif_level,
    break_in_plot = cfg$plots$waterfall$break_in_plot,
    top_padding = cfg$plots$waterfall$top_padding,
    custom_title = cfg$plots$waterfall$custom_title
  )
  return(p)
}