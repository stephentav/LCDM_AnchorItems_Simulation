#!/usr/bin/env Rscript

## Run every validated flexMIRT control file listed in flexmirt_manifest.csv
## through irtQ::run_flexmirt(). Results are checkpointed after every run.

options(stringsAsFactors = FALSE)

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1L]))) else getwd()
}

SCRIPT_DIR <- script_path()

parse_args <- function(args) {
  config <- list(
    data_dir = file.path(SCRIPT_DIR, "generated_data"),
    flex_dir = NA_character_,
    results_file = "flexmirt_run_results.csv",
    limit = NA_integer_,
    data_condition = NA_character_,
    analysis_condition = NA_character_,
    rep = NA_integer_,
    show_output = FALSE,
    overwrite = FALSE,
    dry_run = FALSE
  )
  for (arg in args) {
    if (arg == "--show-output") config$show_output <- TRUE
    else if (arg == "--overwrite") config$overwrite <- TRUE
    else if (arg == "--dry-run") config$dry_run <- TRUE
    else if (grepl("^--data-dir=", arg)) config$data_dir <- sub("^--data-dir=", "", arg)
    else if (grepl("^--flex-dir=", arg)) config$flex_dir <- sub("^--flex-dir=", "", arg)
    else if (grepl("^--results-file=", arg)) config$results_file <- sub("^--results-file=", "", arg)
    else if (grepl("^--limit=", arg)) config$limit <- as.integer(sub("^--limit=", "", arg))
    else if (grepl("^--data-condition=", arg)) config$data_condition <- sub("^--data-condition=", "", arg)
    else if (grepl("^--analysis-condition=", arg)) config$analysis_condition <- sub("^--analysis-condition=", "", arg)
    else if (grepl("^--rep=", arg)) config$rep <- as.integer(sub("^--rep=", "", arg))
    else if (arg %in% c("--help", "-h")) {
      cat(paste(
        "Usage: Rscript RunFlexmirt.R [options]",
        "  --data-dir=PATH              Generated-data directory",
        "  --flex-dir=PATH              flexMIRT installation directory",
        "  --results-file=NAME          Results CSV within data-dir",
        "  --limit=N                    Run only the first N selected files",
        "  --data-condition=ID          Restrict to one data condition",
        "  --analysis-condition=ID      Restrict to one analysis condition",
        "  --rep=N                      Restrict to one replication",
        "  --show-output                Show flexMIRT console output",
        "  --overwrite                  Rerun models already recorded",
        "  --dry-run                    Validate inputs without launching flexMIRT",
        sep = "\n"
      ), "\n")
      quit(save = "no", status = 0L)
    } else stop("Unknown argument: ", arg)
  }
  if (!is.na(config$limit) && config$limit < 1L) stop("--limit must be >= 1")
  if (!is.na(config$rep) && config$rep < 1L) stop("--rep must be >= 1")
  config$data_dir <- normalizePath(config$data_dir, mustWork = FALSE)
  if (!is.na(config$flex_dir)) config$flex_dir <- normalizePath(config$flex_dir, mustWork = FALSE)
  config
}

read_csv <- function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, quote = TRUE, na = "NA")

first_match <- function(pattern, lines, replacement = "\\1", default = NA_character_) {
  match <- grep(pattern, lines, value = TRUE, ignore.case = TRUE, perl = TRUE)
  if (!length(match)) return(default)
  sub(pattern, replacement, match[1L], ignore.case = TRUE, perl = TRUE)
}

parse_flexmirt_output <- function(model_path) {
  base <- sub("\\.flexmirt$", "", model_path, ignore.case = TRUE)
  irt_path <- paste0(base, "-irt.txt")
  exists <- file.exists(irt_path) && !is.na(file.info(irt_path)$size) && file.info(irt_path)$size > 0
  if (!exists) {
    return(list(
      irt_relative_path = NA_character_, output_exists = FALSE,
      normal_termination = FALSE, first_order_converged = FALSE,
      second_order_maximum = NA, cycles = NA_integer_, max_parameter_change = NA_real_,
      information_condition_number = NA_real_, flexmirt_error = "No nonempty -irt.txt output"
    ))
  }
  lines <- readLines(irt_path, warn = FALSE)
  normal <- any(grepl("flexMIRT\\(R\\) engine status:\\s*Normal termination", lines, ignore.case = TRUE))
  first_order <- any(grepl("First-order test:\\s*Convergence criteria satisfied", lines, ignore.case = TRUE))
  not_maximum <- any(grepl("Second-order test:\\s*Solution is not a maximum", lines, ignore.case = TRUE))
  is_maximum <- any(grepl("Second-order test:\\s*Solution is a maximum", lines, ignore.case = TRUE))
  cycles <- suppressWarnings(as.integer(first_match(
    ".*Number of cycles completed:\\s*([0-9]+).*", lines
  )))
  max_change <- suppressWarnings(as.numeric(first_match(
    ".*Maximum parameter change \\(P#\\):\\s*([-+0-9.eE]+).*", lines
  )))
  condition_number <- suppressWarnings(as.numeric(first_match(
    ".*Condition number of information matrix:\\s*([-+0-9.eE]+).*", lines
  )))
  error_lines <- grep("^\\s*(ERROR|FATAL)(:|\\s)", lines, value = TRUE, ignore.case = TRUE)
  list(
    irt_path = irt_path,
    output_exists = TRUE,
    normal_termination = normal,
    first_order_converged = first_order,
    second_order_maximum = if (not_maximum) FALSE else if (is_maximum) TRUE else NA,
    cycles = cycles,
    max_parameter_change = max_change,
    information_condition_number = condition_number,
    flexmirt_error = paste(unique(trimws(error_lines)), collapse = " | ")
  )
}

run_one <- function(task, config) {
  model_path <- normalizePath(file.path(config$data_dir, task$model_relative_path), mustWork = FALSE)
  if (!file.exists(model_path) || basename(model_path) %in% c("", ".", "..") || startsWith(basename(model_path), "._")) {
    stop("Validated model file not found: ", model_path)
  }
  warnings <- character(0)
  execution_error <- ""
  start <- proc.time()[["elapsed"]]
  if (!config$dry_run) {
    old_dir <- getwd()
    on.exit(setwd(old_dir), add = TRUE)
    setwd(dirname(model_path))
    tryCatch(
      withCallingHandlers(
        irtQ::run_flexmirt(
          file.syntax = basename(model_path),
          dir.flex = if (is.na(config$flex_dir)) NULL else config$flex_dir,
          show.output.on.console = config$show_output
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) execution_error <<- conditionMessage(e)
    )
  }
  elapsed <- proc.time()[["elapsed"]] - start
  parsed <- parse_flexmirt_output(model_path)
  irt_relative <- if (!is.null(parsed$irt_path)) {
    substring(parsed$irt_path, nchar(config$data_dir) + 2L)
  } else NA_character_

  data.frame(
    analysis_condition_id = task$analysis_condition_id,
    data_condition_id = task$data_condition_id,
    rep = task$rep,
    Qspec = task$Qspec,
    calibration = task$calibration,
    model_relative_path = task$model_relative_path,
    irt_relative_path = irt_relative,
    execution_attempted = !config$dry_run,
    output_exists = parsed$output_exists,
    normal_termination = parsed$normal_termination,
    first_order_converged = parsed$first_order_converged,
    second_order_maximum = parsed$second_order_maximum,
    cycles = parsed$cycles,
    max_parameter_change = parsed$max_parameter_change,
    information_condition_number = parsed$information_condition_number,
    elapsed_seconds = elapsed,
    runner_warning = paste(unique(warnings), collapse = " | "),
    execution_error = execution_error,
    flexmirt_error = parsed$flexmirt_error,
    irtQ_version = as.character(utils::packageVersion("irtQ")),
    run_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
}

select_tasks <- function(manifest, config) {
  tasks <- manifest
  if (!is.na(config$data_condition)) tasks <- tasks[tasks$data_condition_id == config$data_condition, ]
  if (!is.na(config$analysis_condition)) tasks <- tasks[tasks$analysis_condition_id == config$analysis_condition, ]
  if (!is.na(config$rep)) tasks <- tasks[tasks$rep == config$rep, ]
  tasks <- tasks[order(tasks$analysis_condition_id, tasks$rep), ]
  if (!is.na(config$limit)) tasks <- head(tasks, config$limit)
  rownames(tasks) <- NULL
  tasks
}

main <- function() {
  config <- parse_args(commandArgs(trailingOnly = TRUE))
  if (!requireNamespace("irtQ", quietly = TRUE)) stop("Install the irtQ package before running this script")
  manifest_path <- file.path(config$data_dir, "flexmirt_manifest.csv")
  if (!file.exists(manifest_path)) stop("flexmirt_manifest.csv not found in ", config$data_dir)
  manifest <- read_csv(manifest_path)
  if (!all(manifest$validation_passed)) stop("Manifest contains unvalidated model files")
  tasks <- select_tasks(manifest, config)
  if (!nrow(tasks)) stop("No flexMIRT files match the requested filters")

  missing_models <- !file.exists(file.path(config$data_dir, tasks$model_relative_path)) |
    startsWith(basename(tasks$model_relative_path), "._")
  if (any(missing_models)) stop("Missing or invalid model files: ", sum(missing_models))

  results_path <- file.path(config$data_dir, config$results_file)
  existing <- if (file.exists(results_path)) read_csv(results_path) else NULL
  if (!is.null(existing) && !config$overwrite) {
    completed <- paste(existing$analysis_condition_id, existing$rep, sep = "::")
    keys <- paste(tasks$analysis_condition_id, tasks$rep, sep = "::")
    tasks <- tasks[!keys %in% completed, ]
  }
  if (!nrow(tasks)) {
    message("All selected flexMIRT runs are already recorded in ", results_path)
    return(invisible(NULL))
  }
  if (config$overwrite && !is.null(existing)) {
    rerun <- paste(tasks$analysis_condition_id, tasks$rep, sep = "::")
    existing <- existing[!paste(existing$analysis_condition_id, existing$rep, sep = "::") %in% rerun, ]
  }

  new_results <- vector("list", nrow(tasks))
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, , drop = FALSE]
    message(sprintf("flexMIRT run %d/%d: %s %s rep %d", i, nrow(tasks),
                    task$analysis_condition_id, task$data_condition_id, task$rep))
    new_results[[i]] <- run_one(task, config)
    combined <- do.call(rbind, c(if (!is.null(existing)) list(existing) else list(), new_results[seq_len(i)]))
    write_csv(combined, results_path)
  }
  message(if (config$dry_run) "Dry run complete: " else "flexMIRT batch complete: ",
          nrow(tasks), " files; results in ", results_path)
}

if (sys.nframe() == 0L) main()

