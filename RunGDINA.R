#!/usr/bin/env Rscript

## Fit the generated datasets with the GDINA package using the logit G-DINA
## parameterization (LCDM). The same data, Q specifications, and anchor equality
## structure used by the flexMIRT files are retained.

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
    results_file = "gdina_convergence.csv",
    maxitr = 20000L,
    conv_crit = 1e-4,
    conv_type = "ip,mp",
    nstarts = 3L,
    lower_p = 1e-4,
    upper_p = 0.9999,
    analysis_limit = NA_integer_,
    data_condition = NA_character_,
    rep = NA_integer_,
    save_fits = FALSE,
    overwrite = FALSE
  )
  for (arg in args) {
    if (arg == "--save-fits") config$save_fits <- TRUE
    else if (arg == "--overwrite") config$overwrite <- TRUE
    else if (grepl("^--data-dir=", arg)) config$data_dir <- sub("^--data-dir=", "", arg)
    else if (grepl("^--results-file=", arg)) config$results_file <- sub("^--results-file=", "", arg)
    else if (grepl("^--maxitr=", arg)) config$maxitr <- as.integer(sub("^--maxitr=", "", arg))
    else if (grepl("^--conv-crit=", arg)) config$conv_crit <- as.numeric(sub("^--conv-crit=", "", arg))
    else if (grepl("^--conv-type=", arg)) config$conv_type <- sub("^--conv-type=", "", arg)
    else if (grepl("^--nstarts=", arg)) config$nstarts <- as.integer(sub("^--nstarts=", "", arg))
    else if (grepl("^--analysis-limit=", arg)) config$analysis_limit <- as.integer(sub("^--analysis-limit=", "", arg))
    else if (grepl("^--data-condition=", arg)) config$data_condition <- sub("^--data-condition=", "", arg)
    else if (grepl("^--rep=", arg)) config$rep <- as.integer(sub("^--rep=", "", arg))
    else if (arg %in% c("--help", "-h")) {
      cat(paste(
        "Usage: Rscript RunGDINA.R [options]",
        "  --data-dir=PATH          Generated-data directory",
        "  --results-file=NAME      Results CSV within data-dir",
        "  --maxitr=N               Maximum EM cycles (default: 20000)",
        "  --conv-crit=VALUE        Convergence criterion (default: 1e-4)",
        "  --conv-type=TYPES        Comma-separated criteria (default: ip,mp)",
        "  --nstarts=N              Random starting sets (default: 3)",
        "  --analysis-limit=N       Run only the first N selected fits",
        "  --data-condition=ID      Restrict to one data condition",
        "  --rep=N                  Restrict to one replication number",
        "  --save-fits              Save fitted GDINA objects under gdina_fits/",
        "  --overwrite              Rerun fits already present in results CSV",
        sep = "\n"
      ), "\n")
      quit(save = "no", status = 0L)
    } else stop("Unknown argument: ", arg)
  }
  if (is.na(config$maxitr) || config$maxitr < 1L) stop("--maxitr must be >= 1")
  if (is.na(config$nstarts) || config$nstarts < 1L) stop("--nstarts must be >= 1")
  if (!is.finite(config$conv_crit) || config$conv_crit <= 0) stop("--conv-crit must be positive")
  config$conv_type <- trimws(strsplit(config$conv_type, ",", fixed = TRUE)[[1L]])
  allowed <- c("ip", "mp", "delta", "neg2LL", "neg2LLratio")
  if (!all(config$conv_type %in% allowed)) stop("Unsupported --conv-type value")
  config$data_dir <- normalizePath(config$data_dir, mustWork = FALSE)
  config
}

read_csv <- function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, quote = TRUE, na = "NA")

read_q <- function(data_dir, qspec) {
  q <- read_csv(file.path(data_dir, paste0("Q_", qspec, ".csv")))
  result <- as.matrix(q[, -1L, drop = FALSE])
  storage.mode(result) <- "integer"
  rownames(result) <- q[[1L]]
  result
}

read_anchor_items <- function(path) {
  x <- read_csv(path)
  if (!nrow(x)) integer(0) else as.integer(x[[1L]])
}

make_union_find <- function(keys) {
  parent <- stats::setNames(keys, keys)
  find_root <- function(x) {
    while (parent[[x]] != x) x <- parent[[x]]
    x
  }
  unite <- function(a, b) {
    root_a <- find_root(a); root_b <- find_root(b)
    if (root_a != root_b) parent[[root_b]] <<- root_a
  }
  roots <- function() vapply(keys, find_root, character(1))
  list(unite = unite, roots = roots)
}

build_measurement_design <- function(rep_dir, q_matrix) {
  grades <- 3:5
  n <- nrow(read_csv(file.path(rep_dir, "responses_G3.csv")))
  responses <- lapply(grades, function(grade) {
    x <- read_csv(file.path(rep_dir, paste0("responses_G", grade, ".csv")))
    if (nrow(x) != n || ncol(x) != nrow(q_matrix)) stop("Unexpected response dimensions in ", rep_dir)
    as.matrix(x)
  })
  keys <- unlist(lapply(grades, function(grade) paste0("G", grade, "_I", seq_len(nrow(q_matrix)))))
  uf <- make_union_find(keys)
  anchors34 <- read_anchor_items(file.path(rep_dir, "anchors_G3_G4.csv"))
  anchors45 <- read_anchor_items(file.path(rep_dir, "anchors_G4_G5.csv"))
  for (item in anchors34) uf$unite(paste0("G3_I", item), paste0("G4_I", item))
  for (item in anchors45) uf$unite(paste0("G4_I", item), paste0("G5_I", item))
  roots <- uf$roots()
  unique_roots <- unique(roots)
  root_column <- stats::setNames(seq_along(unique_roots), unique_roots)
  key_column <- unname(root_column[roots])

  dat <- matrix(NA_integer_, nrow = 3L * n, ncol = length(unique_roots))
  q_out <- matrix(NA_integer_, nrow = length(unique_roots), ncol = ncol(q_matrix))
  source_item <- integer(length(unique_roots))
  for (grade_index in seq_along(grades)) {
    rows <- (grade_index - 1L) * n + seq_len(n)
    for (item in seq_len(nrow(q_matrix))) {
      key_index <- (grade_index - 1L) * nrow(q_matrix) + item
      column <- key_column[key_index]
      dat[rows, column] <- responses[[grade_index]][, item]
      if (all(is.na(q_out[column, ]))) {
        q_out[column, ] <- q_matrix[item, ]
        source_item[column] <- item
      } else if (!identical(as.integer(q_out[column, ]), as.integer(q_matrix[item, ]))) {
        stop("An anchor equality class has inconsistent Q rows")
      }
    }
  }
  colnames(dat) <- paste0("M", seq_len(ncol(dat)))
  rownames(q_out) <- colnames(dat)
  colnames(q_out) <- colnames(q_matrix)
  list(
    dat = dat,
    Q = q_out,
    group = factor(rep(paste0("G", grades), each = n), levels = paste0("G", grades)),
    anchors34 = anchors34,
    anchors45 = anchors45,
    n_per_grade = n,
    n_measurement_parameters = ncol(dat)
  )
}

safe_extract <- function(fit, what, default = NA) {
  tryCatch(GDINA::extract(fit, what), error = function(e) default)
}

fit_one <- function(task, config) {
  rep_dir <- file.path(config$data_dir, task$relative_path)
  q_matrix <- read_q(config$data_dir, task$Qspec)
  design <- build_measurement_design(rep_dir, q_matrix)
  group_arg <- if (task$calibration == "multigroup") design$group else NULL
  random_seed <- 700000L + as.integer(sub("A", "", task$analysis_condition_id, fixed = TRUE)) * 100L + task$rep
  warnings <- character(0)
  start_time <- proc.time()[["elapsed"]]

  fit <- tryCatch(
    withCallingHandlers(
      GDINA::GDINA(
        dat = design$dat,
        Q = design$Q,
        model = "logitGDINA",
        att.dist = "saturated",
        mono.constraint = FALSE,
        group = group_arg,
        verbose = 0,
        control = list(
          maxitr = config$maxitr,
          conv.crit = config$conv_crit,
          conv.type = config$conv_type,
          nstarts = config$nstarts,
          lower.p = config$lower_p,
          upper.p = config$upper_p,
          randomseed = random_seed,
          MstepMessage = 0
        )
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[["elapsed"]] - start_time

  if (inherits(fit, "error")) {
    return(data.frame(
      analysis_condition_id = task$analysis_condition_id,
      data_condition_id = task$data_condition_id,
      rep = task$rep, Qspec = task$Qspec, calibration = task$calibration,
      N_per_grade = task$N_per_grade, N_total = 3L * task$N_per_grade,
      n_groups = if (task$calibration == "multigroup") 3L else 1L,
      n_measurement_parameters = design$n_measurement_parameters,
      n_anchors_G3_G4 = length(design$anchors34), n_anchors_G4_G5 = length(design$anchors45),
      converged = FALSE, iterations = NA_integer_, max_iterations_reached = NA,
      logLik = NA_real_, deviance = NA_real_, AIC = NA_real_, BIC = NA_real_,
      elapsed_seconds = elapsed, random_seed = random_seed,
      warning = paste(unique(warnings), collapse = " | "), error = conditionMessage(fit),
      GDINA_version = as.character(utils::packageVersion("GDINA")),
      stringsAsFactors = FALSE
    ))
  }

  converged <- isTRUE(safe_extract(fit, "convergence", FALSE))
  iterations <- as.integer(safe_extract(fit, "nitr", NA_integer_))[1L]
  result <- data.frame(
    analysis_condition_id = task$analysis_condition_id,
    data_condition_id = task$data_condition_id,
    rep = task$rep, Qspec = task$Qspec, calibration = task$calibration,
    N_per_grade = task$N_per_grade, N_total = 3L * task$N_per_grade,
    n_groups = if (task$calibration == "multigroup") 3L else 1L,
    n_measurement_parameters = design$n_measurement_parameters,
    n_anchors_G3_G4 = length(design$anchors34), n_anchors_G4_G5 = length(design$anchors45),
    converged = converged, iterations = iterations,
    max_iterations_reached = !converged && !is.na(iterations) && iterations >= config$maxitr,
    logLik = as.numeric(stats::logLik(fit)), deviance = as.numeric(stats::deviance(fit)),
    AIC = as.numeric(stats::AIC(fit)), BIC = as.numeric(stats::BIC(fit)),
    elapsed_seconds = elapsed, random_seed = random_seed,
    warning = paste(unique(warnings), collapse = " | "), error = "",
    GDINA_version = as.character(utils::packageVersion("GDINA")),
    stringsAsFactors = FALSE
  )
  if (config$save_fits) {
    fit_dir <- file.path(config$data_dir, "gdina_fits", task$data_condition_id, sprintf("rep_%03d", task$rep))
    dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(fit, file.path(fit_dir, paste0("fit_", task$analysis_condition_id, ".rds")))
  }
  result
}

build_tasks <- function(config) {
  analyses <- read_csv(file.path(config$data_dir, "analysis_conditions.csv"))
  reps <- read_csv(file.path(config$data_dir, "rep_manifest_all.csv"))
  tasks <- merge(reps, analyses[, c("analysis_condition_id", "data_condition_id", "Qspec", "calibration")],
                 by = "data_condition_id", sort = FALSE)
  tasks <- tasks[order(tasks$analysis_condition_id, tasks$rep), ]
  if (!is.na(config$data_condition)) tasks <- tasks[tasks$data_condition_id == config$data_condition, ]
  if (!is.na(config$rep)) tasks <- tasks[tasks$rep == config$rep, ]
  if (!is.na(config$analysis_limit)) tasks <- head(tasks, config$analysis_limit)
  rownames(tasks) <- NULL
  tasks
}

main <- function() {
  if (!requireNamespace("GDINA", quietly = TRUE)) stop("Install the GDINA package before running this script")
  config <- parse_args(commandArgs(trailingOnly = TRUE))
  tasks <- build_tasks(config)
  if (!nrow(tasks)) stop("No fits match the requested filters")
  results_path <- file.path(config$data_dir, config$results_file)
  existing <- if (file.exists(results_path)) read_csv(results_path) else NULL
  if (!is.null(existing) && !config$overwrite) {
    completed_keys <- paste(existing$analysis_condition_id, existing$rep, sep = "::")
    task_keys <- paste(tasks$analysis_condition_id, tasks$rep, sep = "::")
    tasks <- tasks[!task_keys %in% completed_keys, ]
  }
  if (!nrow(tasks)) {
    message("All selected GDINA fits are already recorded in ", results_path)
    return(invisible(NULL))
  }

  if (config$overwrite && !is.null(existing)) {
    rerun_keys <- paste(tasks$analysis_condition_id, tasks$rep, sep = "::")
    existing <- existing[!paste(existing$analysis_condition_id, existing$rep, sep = "::") %in% rerun_keys, ]
  }

  new_results <- vector("list", nrow(tasks))
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, , drop = FALSE]
    message(sprintf("GDINA fit %d/%d: %s %s rep %d", i, nrow(tasks),
                    task$analysis_condition_id, task$data_condition_id, task$rep))
    new_results[[i]] <- fit_one(task, config)
    combined <- do.call(rbind, c(if (!is.null(existing)) list(existing) else list(), new_results[seq_len(i)]))
    write_csv(combined, results_path)
  }
  message("GDINA fitting complete: ", nrow(tasks), " new fits; results in ", results_path)
}

if (sys.nframe() == 0L) main()

