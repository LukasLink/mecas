# R/R_snake_functions/bootstrap.R

bootstrap_pipeline <- function(project_root_dir) {
  stopifnot(is.character(project_root_dir), length(project_root_dir) == 1)
  
  # Global R options
  options(bitmapType = "cairo")
  
  if (requireNamespace("knitr", quietly = TRUE)) {
    knitr::opts_chunk$set(echo = FALSE)
  }
  
  cran_packages <- c(
    "tidyverse",
    "Matrix",
    "conflicted",
    "ggplot2",
    "ggrepel",
    "writexl",
    "stringr",
    "ggbreak",
    "yaml",
    "tools",
    "logger"
  )
  
  bioc_packages <- c(
    "AnnotationDbi",
    "org.Hs.eg.db"
  )
  
  other_packages <- c(
    "MAUDE"
  )
  
  check_packages <- function(packages, source = "CRAN") {
    missing <- packages[
      !vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
    ]
    
    if (length(missing) > 0) {
      message("Missing ", source, " package(s): ", paste(missing, collapse = ", "))
      
      if (source == "CRAN") {
        message("Install with:")
        message("install.packages(c(", paste(sprintf('"%s"', missing), collapse = ", "), "))")
      }
      
      if (source == "Bioconductor") {
        message("Install with:")
        message("BiocManager::install(c(", paste(sprintf('"%s"', missing), collapse = ", "), "))")
      }
      
      stop("Required package(s) missing.", call. = FALSE)
    }
  }
  
  check_packages(cran_packages, "CRAN")
  check_packages(bioc_packages, "Bioconductor")
  check_packages(other_packages, "other source")
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(Matrix)
    library(conflicted)
    library(MAUDE)
    library(ggplot2)
    library(ggrepel)
    library(writexl)
    library(stringr)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(ggbreak)
    library(yaml)
    library(tools)
    library(logger)
  })
  
  suppressMessages({
    conflicted::conflicts_prefer(dplyr::rename)
    conflicted::conflicts_prefer(dplyr::filter)
    conflicted::conflicts_prefer(dplyr::select)
    conflicted::conflicts_prefer(dplyr::slice)
    conflicted::conflicts_prefer(dplyr::first)
    conflicted::conflicts_prefer(dplyr::desc)
    conflicted::conflicts_prefer(base::setdiff)
    conflicted::conflicts_prefer(base::intersect)
    conflicted::conflicts_prefer(base::unname)
    conflicted::conflicts_prefer(base::setequal)
  })
  
  # Source all pipeline functions

  source(file.path(project_root_dir, "R", "R_snake_functions", "zzz_snake_source_all.R"))
  source_snake_functions(project_root_dir)
  
  invisible(TRUE)
}