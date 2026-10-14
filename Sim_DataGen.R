#!/usr/bin/env Rscript

## Reproducible data generator for the LCDM vertical-linking simulation.
##
## Data are generated once for each combination of sample size, cross-grade
## mastery-prevalence pattern, anchor design, and difficulty pattern. The
## estimation factors (Q specification and calibration method) are recorded in
## a separate analysis manifest so that all four analyses use identical data.

## Test with 5 reps in each condition first

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  config <- list(
    out_dir = file.path(getwd(), "generated_data"),
    n_reps = 5,
    workers = 1L,
    condition_limit = NA_integer_,
    master_seed = 1234L,
    overwrite = FALSE
  )

  for (arg in args) {
    if (arg == "--overwrite") {
      config$overwrite <- TRUE
    } else if (grepl("^--out-dir=", arg)) {
      config$out_dir <- sub("^--out-dir=", "", arg)
    } else if (grepl("^--n-reps=", arg)) {
      config$n_reps <- as.integer(sub("^--n-reps=", "", arg))
    } else if (grepl("^--workers=", arg)) {
      config$workers <- as.integer(sub("^--workers=", "", arg))
    } else if (grepl("^--condition-limit=", arg)) {
      config$condition_limit <- as.integer(sub("^--condition-limit=", "", arg))
    } else if (grepl("^--master-seed=", arg)) {
      config$master_seed <- as.integer(sub("^--master-seed=", "", arg))
    } else if (arg %in% c("--help", "-h")) {
      cat(paste(
        "Usage: Rscript Sim_DataGen.R [options]",
        "  --out-dir=PATH          Output directory (default: ./generated_data)",
        "  --n-reps=N              Replications per data condition (default: 100)",
        "  --workers=N             Parallel workers (default: 1)",
        "  --condition-limit=N     Run only the first N conditions (testing)",
        "  --master-seed=N         Master seed recorded in manifests (default: 1234)",
        "  --overwrite             Replace existing replication directories",
        sep = "\n"
      ), "\n")
      quit(save = "no", status = 0L)
    } else {
      stop("Unknown argument: ", arg)
    }
  }

  if (is.na(config$n_reps) || config$n_reps < 1L) stop("--n-reps must be >= 1")
  if (is.na(config$workers) || config$workers < 1L) stop("--workers must be >= 1")
  if (is.na(config$master_seed) || config$master_seed < 1L) stop("--master-seed must be >= 1")
  config$out_dir <- normalizePath(config$out_dir, mustWork = FALSE)
  config
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE, na = "NA")
}

invlogit <- function(x) stats::plogis(x)

Q_TRUE <- matrix(c(
  1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1,
  1,1,0,0, 1,0,1,0, 1,0,0,1, 0,1,0,1,
  0,1,1,0, 0,0,1,1, 1,0,0,0, 0,1,0,0,
  0,0,1,0, 0,0,0,1, 1,1,0,0, 1,0,1,0,
  1,0,0,1, 0,1,0,1, 0,1,1,0, 0,0,1,1,
  1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1,
  1,1,0,0, 1,0,1,0, 1,0,0,1, 0,1,0,1,
  0,1,1,0, 0,0,1,1
), nrow = 30L, byrow = TRUE,
dimnames = list(paste0("I", seq_len(30L)), paste0("A", seq_len(4L))))

Q_MISSPEC <- Q_TRUE
Q_MISSPEC[c(5L, 6L, 7L, 8L), ] <- rbind(
  c(0,1,0,0), c(0,0,1,0), c(0,0,0,1), c(0,0,0,1)
)

MASTERY_PROBS <- list(
  no = list(`3` = rep(0.5, 4L), `4` = rep(0.5, 4L), `5` = rep(0.5, 4L)),
  slow = list(
    `3` = c(0.30, 0.20, 0.10, 0.15),
    `4` = c(0.40, 0.30, 0.15, 0.25),
    `5` = c(0.50, 0.40, 0.20, 0.35)
  ),
  high = list(
    `3` = c(0.30, 0.20, 0.10, 0.15),
    `4` = c(0.50, 0.40, 0.25, 0.35),
    `5` = c(0.70, 0.60, 0.40, 0.60)
  )
)

mastery_probs <- function(pattern, grade) {
  probs <- MASTERY_PROBS[[pattern]][[as.character(grade)]]
  if (is.null(probs)) stop("Unknown mastery pattern/grade: ", pattern, "/", grade)
  probs
}

simulate_profiles <- function(n, probs) {
  result <- vapply(probs, function(p) stats::rbinom(n, 1L, p), integer(n))
  colnames(result) <- colnames(Q_TRUE)
  result
}

generate_item_params <- function(q_matrix, grade, difficulty_pattern) {
  mu <- if (difficulty_pattern == "constant") {
    -1.5
  } else {
    unname(c(`3` = -1.0, `4` = -1.5, `5` = -2.0)[as.character(grade)])
  }
  if (is.na(mu)) stop("Unsupported grade: ", grade)

  data.frame(
    item = seq_len(nrow(q_matrix)),
    grade = grade,
    delta01 = pmin(stats::rnorm(nrow(q_matrix), mu, 0.25), -0.2),
    main = stats::runif(nrow(q_matrix), 1.0, 2.0),
    twoway = pmax(stats::rnorm(nrow(q_matrix), 0.5, 0.2), 0),
    stringsAsFactors = FALSE
  )
}

item_complexity <- function(q_matrix) rowSums(q_matrix)

choose_balanced <- function(candidates, q_matrix, difficulty, n_pick, initial_coverage = NULL) {
  if (n_pick == 0L) return(integer(0))
  if (length(candidates) < n_pick) stop("Not enough candidate anchor items")
  coverage <- if (is.null(initial_coverage)) rep(0L, ncol(q_matrix)) else initial_coverage
  selected <- integer(0)

  while (length(selected) < n_pick) {
    available <- setdiff(candidates, selected)
    scores <- vapply(available, function(item) {
      required <- which(q_matrix[item, ] == 1L)
      sum(1 / (coverage[required] + 1))
    }, numeric(1))
    best_score <- max(scores)
    tied <- available[scores == best_score]
    chosen <- tied[which.min(difficulty[tied])]
    selected <- c(selected, chosen)
    coverage <- coverage + q_matrix[chosen, ]
  }
  selected
}

select_anchor_items <- function(params, q_matrix, n_anchor, location_mode) {
  n_anchor <- as.integer(n_anchor)
  if (n_anchor == 0L) return(integer(0))
  if (!location_mode %in% c("single", "mixed")) stop("Unknown anchor location: ", location_mode)

  complexity <- item_complexity(q_matrix)
  singles <- which(complexity == 1L)
  doubles <- which(complexity == 2L)

  if (location_mode == "single") {
    selected <- choose_balanced(singles, q_matrix, params$delta01, n_anchor)
  } else {
    n_single <- n_anchor %/% 2L
    n_double <- n_anchor - n_single
    selected_single <- choose_balanced(singles, q_matrix, params$delta01, n_single)
    coverage <- if (length(selected_single)) colSums(q_matrix[selected_single, , drop = FALSE]) else rep(0L, ncol(q_matrix))
    selected_double <- choose_balanced(doubles, q_matrix, params$delta01, n_double, coverage)
    selected <- c(selected_single, selected_double)
  }

  sort(as.integer(selected))
}

generate_multigrade_params <- function(q_matrix, anchors_per_attribute,
                                       anchor_location, difficulty_pattern) {
  n_anchor <- as.integer(anchors_per_attribute) * ncol(q_matrix)

  grade3 <- generate_item_params(q_matrix, 3L, difficulty_pattern)
  anchors34 <- select_anchor_items(grade3, q_matrix, n_anchor, anchor_location)

  grade4 <- generate_item_params(q_matrix, 4L, difficulty_pattern)
  if (length(anchors34)) {
    grade4[anchors34, c("delta01", "main", "twoway")] <-
      grade3[anchors34, c("delta01", "main", "twoway")]
  }
  anchors45 <- select_anchor_items(grade4, q_matrix, n_anchor, anchor_location)

  grade5 <- generate_item_params(q_matrix, 5L, difficulty_pattern)
  if (length(anchors45)) {
    grade5[anchors45, c("delta01", "main", "twoway")] <-
      grade4[anchors45, c("delta01", "main", "twoway")]
  }

  grade3$anchor34 <- grade3$item %in% anchors34
  grade3$anchor45 <- FALSE
  grade4$anchor34 <- grade4$item %in% anchors34
  grade4$anchor45 <- grade4$item %in% anchors45
  grade5$anchor34 <- FALSE
  grade5$anchor45 <- grade5$item %in% anchors45

  list(
    grade3 = grade3,
    grade4 = grade4,
    grade5 = grade5,
    anchors34 = anchors34,
    anchors45 = anchors45
  )
}

simulate_lcdm_responses <- function(alpha, q_matrix, params) {
  n <- nrow(alpha)
  responses <- matrix(0L, nrow = n, ncol = nrow(q_matrix))

  for (item in seq_len(nrow(q_matrix))) {
    required <- which(q_matrix[item, ] == 1L)
    main_term <- if (length(required)) rowSums(alpha[, required, drop = FALSE]) else rep(0, n)
    interaction_term <- rep(0, n)
    if (length(required) >= 2L) {
      pairs <- utils::combn(required, 2L)
      interaction_term <- rowSums(
        alpha[, pairs[1L, ], drop = FALSE] * alpha[, pairs[2L, ], drop = FALSE]
      )
    }
    eta <- params$delta01[item] + params$main[item] * main_term +
      params$twoway[item] * interaction_term
    responses[, item] <- stats::rbinom(n, 1L, invlogit(eta))
  }
  colnames(responses) <- rownames(q_matrix)
  responses
}

build_data_conditions <- function() {
  conditions <- expand.grid(
    N = c(500L, 1000L, 2000L),
    prevalence_pattern = c("no", "slow", "high"),
    anchors_per_attribute = c(0L, 1L, 2L),
    anchor_location = c("single", "mixed"),
    difficulty_pattern = c("constant", "progressive"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  # Anchor placement is undefined when no anchors are present. Retain one
  # zero-anchor row and label its placement explicitly rather than duplicating
  # otherwise identical "single" and "mixed" conditions.
  conditions <- conditions[
    conditions$anchors_per_attribute > 0L |
      conditions$anchor_location == "single",
  ]
  conditions$anchor_location[conditions$anchors_per_attribute == 0L] <- "none"
  conditions <- conditions[order(
    conditions$N, conditions$prevalence_pattern,
    conditions$anchors_per_attribute, conditions$anchor_location,
    conditions$difficulty_pattern
  ), ]
  rownames(conditions) <- NULL
  conditions$data_condition_id <- sprintf("D%04d", seq_len(nrow(conditions)))
  conditions[, c("data_condition_id", setdiff(names(conditions), "data_condition_id"))]
}

build_analysis_conditions <- function(data_conditions) {
  analysis <- merge(
    data_conditions,
    expand.grid(
      Qspec = c("true", "misspec"),
      calibration = c("multigroup", "pooled"),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ),
    by = NULL
  )
  analysis <- analysis[order(
    analysis$data_condition_id, analysis$Qspec, analysis$calibration
  ), ]
  rownames(analysis) <- NULL
  analysis$analysis_condition_id <- sprintf("A%04d", seq_len(nrow(analysis)))
  analysis[, c("analysis_condition_id", "data_condition_id", "Qspec", "calibration",
               setdiff(names(analysis), c("analysis_condition_id", "data_condition_id", "Qspec", "calibration")))]
}

replication_seed <- function(master_seed, condition_index, rep) {
  seed <- as.double(master_seed) + as.double(condition_index) * 100000 + as.double(rep)
  if (seed > .Machine$integer.max) stop("Derived seed exceeds R integer range")
  as.integer(seed)
}

validate_replication <- function(rep_dir, condition, params) {
  expected_anchors <- as.integer(condition$anchors_per_attribute) * ncol(Q_TRUE)
  stopifnot(
    length(params$anchors34) == expected_anchors,
    length(params$anchors45) == expected_anchors,
    !anyDuplicated(params$anchors34),
    !anyDuplicated(params$anchors45)
  )

  coefficient_names <- c("delta01", "main", "twoway")
  if (length(params$anchors34)) {
    stopifnot(isTRUE(all.equal(
      params$grade3[params$anchors34, coefficient_names],
      params$grade4[params$anchors34, coefficient_names],
      check.attributes = FALSE
    )))
  }
  if (length(params$anchors45)) {
    stopifnot(isTRUE(all.equal(
      params$grade4[params$anchors45, coefficient_names],
      params$grade5[params$anchors45, coefficient_names],
      check.attributes = FALSE
    )))
  }

  required_files <- c(
    "anchors_G3_G4.csv", "anchors_G4_G5.csv",
    paste0("item_params_G", 3:5, ".csv"),
    paste0("alpha_G", 3:5, ".csv"),
    paste0("responses_G", 3:5, ".csv")
  )
  stopifnot(all(file.exists(file.path(rep_dir, required_files))))
  TRUE
}

generate_replication <- function(condition, condition_index, rep, config) {
  seed <- replication_seed(config$master_seed, condition_index, rep)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

  condition_dir <- file.path(config$out_dir, condition$data_condition_id)
  rep_dir <- file.path(condition_dir, sprintf("rep_%03d", rep))
  if (dir.exists(rep_dir)) {
    if (!config$overwrite) stop("Replication directory already exists: ", rep_dir)
    unlink(rep_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

  params <- generate_multigrade_params(
    Q_TRUE,
    condition$anchors_per_attribute,
    condition$anchor_location,
    condition$difficulty_pattern
  )

  write_csv(data.frame(item = params$anchors34), file.path(rep_dir, "anchors_G3_G4.csv"))
  write_csv(data.frame(item = params$anchors45), file.path(rep_dir, "anchors_G4_G5.csv"))

  for (grade in 3:5) {
    grade_key <- paste0("grade", grade)
    probs <- mastery_probs(condition$prevalence_pattern, grade)
    alpha <- simulate_profiles(condition$N, probs)
    responses <- simulate_lcdm_responses(alpha, Q_TRUE, params[[grade_key]])

    write_csv(params[[grade_key]], file.path(rep_dir, paste0("item_params_G", grade, ".csv")))
    write_csv(as.data.frame(alpha), file.path(rep_dir, paste0("alpha_G", grade, ".csv")))
    write_csv(as.data.frame(responses), file.path(rep_dir, paste0("responses_G", grade, ".csv")))
  }

  validate_replication(rep_dir, condition, params)

  data.frame(
    data_condition_id = condition$data_condition_id,
    rep = rep,
    seed = seed,
    N_per_grade = condition$N,
    prevalence_pattern = condition$prevalence_pattern,
    anchors_per_attribute = condition$anchors_per_attribute,
    n_anchors_G3_G4 = length(params$anchors34),
    n_anchors_G4_G5 = length(params$anchors45),
    anchor_location = condition$anchor_location,
    difficulty_pattern = condition$difficulty_pattern,
    relative_path = file.path(condition$data_condition_id, sprintf("rep_%03d", rep)),
    validation_passed = TRUE,
    stringsAsFactors = FALSE
  )
}

run_condition <- function(condition_index, data_conditions, config) {
  condition <- data_conditions[condition_index, , drop = FALSE]
  condition_dir <- file.path(config$out_dir, condition$data_condition_id)
  dir.create(condition_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- lapply(seq_len(config$n_reps), function(rep) {
    generate_replication(condition, condition_index, rep, config)
  })
  manifest <- do.call(rbind, rows)
  write_csv(manifest, file.path(condition_dir, "rep_manifest.csv"))
  manifest
}

main <- function() {
  config <- parse_args(commandArgs(trailingOnly = TRUE))
  if (dir.exists(config$out_dir) && length(list.files(config$out_dir, all.files = TRUE, no.. = TRUE)) && !config$overwrite) {
    stop("Output directory is not empty. Use a new --out-dir or pass --overwrite: ", config$out_dir)
  }
  dir.create(config$out_dir, recursive = TRUE, showWarnings = FALSE)

  data_conditions <- build_data_conditions()
  analysis_conditions <- build_analysis_conditions(data_conditions)
  if (!is.na(config$condition_limit)) {
    if (config$condition_limit < 1L) stop("--condition-limit must be >= 1")
    data_conditions_to_run <- head(data_conditions, config$condition_limit)
  } else {
    data_conditions_to_run <- data_conditions
  }

  write_csv(data_conditions, file.path(config$out_dir, "data_conditions.csv"))
  write_csv(analysis_conditions, file.path(config$out_dir, "analysis_conditions.csv"))
  write_csv(cbind(item = rownames(Q_TRUE), as.data.frame(Q_TRUE)), file.path(config$out_dir, "Q_true.csv"))
  write_csv(cbind(item = rownames(Q_MISSPEC), as.data.frame(Q_MISSPEC)), file.path(config$out_dir, "Q_misspec.csv"))
  write_csv(data.frame(
    setting = c("master_seed", "n_reps", "workers", "condition_limit", "rng_kind", "normal_kind", "sample_kind"),
    value = c(config$master_seed, config$n_reps, config$workers,
              ifelse(is.na(config$condition_limit), "all", config$condition_limit),
              "Mersenne-Twister", "Inversion", "Rejection")
  ), file.path(config$out_dir, "run_config.csv"))
  capture.output(sessionInfo(), file = file.path(config$out_dir, "sessionInfo.txt"))

  indices <- seq_len(nrow(data_conditions_to_run))
  if (config$workers == 1L) {
    manifests <- lapply(indices, function(index) {
      message(sprintf("Generating %s (%d/%d)", data_conditions_to_run$data_condition_id[index], index, length(indices)))
      run_condition(index, data_conditions_to_run, config)
    })
  } else if (.Platform$OS.type == "unix") {
    manifests <- parallel::mclapply(
      indices,
      function(index) run_condition(index, data_conditions_to_run, config),
      mc.cores = config$workers,
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
  } else {
    warning("Parallel generation is unavailable on this platform; using one worker")
    manifests <- lapply(indices, function(index) {
      run_condition(index, data_conditions_to_run, config)
    })
  }

  all_manifest <- do.call(rbind, manifests)
  write_csv(all_manifest, file.path(config$out_dir, "rep_manifest_all.csv"))
  message("Generation complete: ", nrow(all_manifest), " replications in ", config$out_dir)
}

if (sys.nframe() == 0L) main()
