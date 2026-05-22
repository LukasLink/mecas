# 22_run_shell_step.R
#-------------------------------------------------------------------------------
# Run shell commands locally or through SLURM
#-------------------------------------------------------------------------------

run_shell_step <- function(step_name,
                           script_path,
                           slurm_settings,
                           args = character(),
                           machine = "local",
                           log_dir = NULL,
                           extra_slurm_sbatch_lines = character()
                           ) {
  
  if (!is.character(step_name) || length(step_name) != 1 || !nzchar(step_name)) {
    stop_log("`step_name` must be a single non-empty string.")
  }
  
  if (!file.exists(script_path)) {
    stop_log("Shell script does not exist: ", script_path)
  }
  if (!machine %in% c("local", "slurm")) {
    stop_log(
      "`machine` must be either 'local' or 'slurm'. Current value: ",
      machine
    )
  }
  
  if (is.null(log_dir)) {
    log_dir <- getwd()
  }
  
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (machine == "local") {
    return(run_shell_step_local(
      step_name = step_name,
      script_path = script_path,
      args = args,
      log_dir = log_dir
    ))
  }
  
  if (machine == "slurm") {
    return(run_shell_step_slurm(
      step_name = step_name,
      script_path = script_path,
      args = args,
      log_dir = log_dir,
      slurm_settings = slurm_settings,
      extra_slurm_sbatch_lines = extra_slurm_sbatch_lines
    ))
  }
}
run_shell_step_local <- function(step_name,
                                 script_path,
                                 args = character(),
                                 log_dir) {
  
  stdout_log <- file.path(log_dir, paste0(step_name, ".local.out.log"))
  stderr_log <- file.path(log_dir, paste0(step_name, ".local.err.log"))
  
  args <- as.character(args)
  script_path <- normalizePath(script_path, mustWork = TRUE)
  
  
  logger::log_info("Running local shell step: {step_name}")
  logger::log_info("Script: {script_path}")
  logger::log_info("stdout log: {stdout_log}")
  logger::log_info("stderr log: {stderr_log}")
  
  status <- system2(
    command = "bash",
    args = c(script_path, args),
    stdout = stdout_log,
    stderr = stderr_log
  )
  
  if (!identical(status, 0L)) {
    stop_log(
      "Shell step failed: ", step_name, "\n",
      "Exit status: ", status, "\n",
      "stderr log: ", stderr_log
    )
  }
  
  logger::log_info("Finished local shell step: {step_name}")
  
  invisible(list(
    machine = "local",
    step_name = step_name,
    script_path = script_path,
    args = args,
    status = status,
    stdout_log = stdout_log,
    stderr_log = stderr_log
  ))
}

parse_sbatch_job_id <- function(sbatch_output) {
  output <- paste(sbatch_output, collapse = " ")
  
  match <- regmatches(
    output,
    regexpr("[0-9]+", output)
  )
  
  if (length(match) == 0 || is.na(match) || !nzchar(match)) {
    warning(
      "Could not parse SLURM job ID from sbatch output: ",
      output,
      call. = FALSE
    )
    return(NA_character_)
  }
  
  match
}

run_shell_step_slurm <- function(step_name,
                                 script_path,
                                 args = character(),
                                 log_dir,
                                 slurm_settings,
                                 extra_slurm_sbatch_lines) {
  
  args <- as.character(args)
  script_path <- normalizePath(script_path, mustWork = TRUE)
  
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  
  account   <- slurm_settings$account
  qos       <- slurm_settings$qos
  cpus      <- slurm_settings$cpus
  mem       <- slurm_settings$mem
  wall_time <- slurm_settings$wall_time
  partition <- slurm_settings$partition
  array     <- slurm_settings$array
  email     <- slurm_settings$email
  
  stdout_log <- file.path(log_dir, paste0(step_name, "_slurm_%x_%A_%a.out.log"))
  stderr_log <- file.path(log_dir, paste0(step_name, "_slurm_%x_%A_%a.err.log"))
  
  sbatch_script <- file.path(log_dir, paste0(step_name, ".sbatch.sh"))
  
  sbatch_lines <- c(
    "#!/usr/bin/env bash",
    paste0("#SBATCH --job-name=", step_name),
    paste0("#SBATCH --cpus-per-task=", cpus),
    paste0("#SBATCH --mem=", mem),
    paste0("#SBATCH --time=", wall_time),
    paste0("#SBATCH --output=", stdout_log),
    paste0("#SBATCH --error=", stderr_log)
  )
  
  if (!is.null(account) && !is.na(account) && nzchar(account)) {
    sbatch_lines <- c(
      sbatch_lines,
      paste0("#SBATCH --account=", account)
    )
  }
  
  if (!is.null(qos) && !is.na(qos) && nzchar(qos)) {
    sbatch_lines <- c(
      sbatch_lines,
      paste0("#SBATCH --qos=", qos)
    )
  }
  
  if (!is.null(partition) && !is.na(partition) && nzchar(partition)) {
    sbatch_lines <- c(
      sbatch_lines,
      paste0("#SBATCH --partition=", partition)
    )
  }
  
  if (!is.null(array) && !is.na(array) && as.integer(array) > 1) {
    sbatch_lines <- c(
      sbatch_lines,
      paste0("#SBATCH --array=0-", as.integer(array) - 1)
    )
  }
  
  if (!is.null(email) && !is.na(email) && nzchar(email)) {
    sbatch_lines <- c(
      sbatch_lines,
      paste0("#SBATCH --mail-user=", email),
      "#SBATCH --mail-type=BEGIN,END,FAIL"
    )
  }
  
  if (!is.character(extra_slurm_sbatch_lines)) {
    stop_log("`extra_slurm_sbatch_lines` must be a character vector.")
  }
  
  if (length(extra_slurm_sbatch_lines) > 0) {
    sbatch_lines <- c(
      sbatch_lines,
      extra_slurm_sbatch_lines
    )
  }
  
  command_line <- paste(
    "bash",
    shQuote(script_path),
    paste(shQuote(args), collapse = " ")
  )
  
  body_lines <- c(
    "",
    "set -euo pipefail",
    "",
    "echo \"Started at: $(date)\"",
    "echo \"Running on host: $(hostname)\"",
    "echo \"SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}\"",
    "echo \"SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-NA}\"",
    paste0("echo ", shQuote(paste0("Command: ", command_line))),
    "",
    command_line,
    "",
    "echo \"Finished at: $(date)\""
  )
  
  writeLines(c(sbatch_lines, body_lines), sbatch_script)
  
  logger::log_info("Submitting SLURM shell step: {step_name}")
  logger::log_info("SBATCH script: {sbatch_script}")
  
  sbatch_output <- system2(
    command = "sbatch",
    args = sbatch_script,
    stdout = TRUE,
    stderr = TRUE
  )
  
  sbatch_status <- attr(sbatch_output, "status")
  
  if (!is.null(sbatch_status) && sbatch_status != 0) {
    stop_log(
      "sbatch submission failed for step: ", step_name, "\n",
      "Output:\n",
      paste(sbatch_output, collapse = "\n")
    )
  }
  
  logger::log_info("sbatch output: {paste(sbatch_output, collapse = ' ')}")
  
  job_id <- parse_sbatch_job_id(sbatch_output)
  
  wait_for_slurm_job(
    job_id = job_id,
    poll_interval_seconds = 30
  )
  
  invisible(list(
    machine = "slurm",
    step_name = step_name,
    script_path = script_path,
    args = args,
    sbatch_script = sbatch_script,
    sbatch_output = sbatch_output,
    job_id = job_id,
    stdout_log = stdout_log,
    stderr_log = stderr_log
  ))
}

wait_for_slurm_job <- function(job_id,
                               poll_interval_seconds = 30,
                               message_interval_minutes = 30) {
  
  if (is.na(job_id) || !nzchar(job_id)) {
    stop_log("Cannot wait for SLURM job because `job_id` is missing.")
  }
  
  logger::log_info("Waiting for SLURM job to finish: {job_id}")
  
  seconds_counter <- 0
  logged_message_counter <- 1
  message_interval_seconds <- message_interval_minutes * 60
  repeat {
    squeue_out <- system2(
      command = "squeue",
      args = c("-j", job_id, "-h"),
      stdout = TRUE,
      stderr = TRUE
    )
    
    squeue_status <- attr(squeue_out, "status")
    
    if (!is.null(squeue_status) && squeue_status != 0) {
      stop_log(
        "Failed to query SLURM job with squeue.\n",
        "Job ID: ", job_id, "\n",
        "squeue output:\n",
        paste(squeue_out, collapse = "\n")
      )
    }
    
    # If squeue returns no rows, the job is no longer running/pending.
    if (length(squeue_out) == 0) {
      break
    }
    
    seconds_counter <- seconds_counter + poll_interval_seconds
    if (seconds_counter >= (logged_message_counter * message_interval_seconds)){
      elapsed_minutes <- round(seconds_counter / 60, 1)
      logger::log_info("SLURM job {job_id} still running or pending. Time elapsed {elapsed_minutes} min")
      logged_message_counter <- logged_message_counter + 1
    }
    Sys.sleep(poll_interval_seconds)
  }
  
  logger::log_info("SLURM job left queue: {job_id}")
  
  check_slurm_job_success(job_id)
  
  invisible(TRUE)
}
check_slurm_job_success <- function(job_id) {
  
  sacct_out <- system2(
    command = "sacct",
    args = c(
      "-j", job_id,
      "--format=JobID,State,ExitCode",
      "--parsable2",
      "--noheader",
      "-X"
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  
  sacct_status <- attr(sacct_out, "status")
  
  if (!is.null(sacct_status) && sacct_status != 0) {
    stop_log(
      "SLURM job finished, but final status could not be checked with sacct.\n",
      "Job ID: ", job_id, "\n",
      "sacct output:\n",
      paste(sacct_out, collapse = "\n")
    )
  }
  
  if (length(sacct_out) == 0) {
    stop_log(
      "SLURM job finished, but sacct returned no status information.\n",
      "Job ID: ", job_id
    )
  }
  
  status_df <- do.call(
    rbind,
    strsplit(sacct_out, "\\|", fixed = FALSE)
  )
  
  status_df <- as.data.frame(status_df, stringsAsFactors = FALSE)
  colnames(status_df) <- c("job_id", "state", "exit_code")
  
  bad_states <- status_df[
    !grepl("^COMPLETED", status_df$state),
    ,
    drop = FALSE
  ]
  
  if (nrow(bad_states) > 0) {
    stop_log(
      "SLURM job did not complete successfully.\n",
      "Job ID: ", job_id, "\n",
      "Observed states:\n",
      paste(
        paste(
          bad_states$job_id,
          bad_states$state,
          bad_states$exit_code,
          sep = "\t"
        ),
        collapse = "\n"
      )
    )
  }
  
  logger::log_info("SLURM job completed successfully: {job_id}")
  
  invisible(TRUE)
}