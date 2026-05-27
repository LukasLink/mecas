# functions/zzz_source_all.R
this_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = TRUE))

r_files <- sort(list.files(this_dir, pattern = "\\.R$", full.names = TRUE))
r_files <- r_files[!grepl("zzz_source_all\\.R$", r_files)]


# For debugging
# invisible(lapply(r_files, function(f) {
#   message("Sourcing: ", f)
#   
#   tryCatch(
#     {
#       source(f, local = .GlobalEnv)
#       message("Successfully sourced: ", f)
#     },
#     error = function(e) {
#       message("ERROR while sourcing: ", f)
#       message("Error message: ", conditionMessage(e))
#       stop(e)
#     }
#   )
# }))


invisible(lapply(r_files, source, local = .GlobalEnv))