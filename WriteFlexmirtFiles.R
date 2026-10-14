#!/usr/bin/env Rscript

## Build flexMIRT control files for every generated replication and analysis.
## The four source templates are read but never modified.

options(stringsAsFactors = FALSE)

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1L]))) else getwd()
}

SCRIPT_DIR <- script_path()
source(file.path(SCRIPT_DIR, "UpdateTol.R"))

parse_args <- function(args) {
  config <- list(
    data_dir = file.path(SCRIPT_DIR, "generated_data"),
    templates_dir = SCRIPT_DIR,
    etol = "1e-4",
    mtol = "1e-7",
    maxe = "20000",
    maxm = "100",
    processors = "1",
    overwrite = FALSE
  )
  for (arg in args) {
    if (arg == "--overwrite") config$overwrite <- TRUE
    else if (grepl("^--data-dir=", arg)) config$data_dir <- sub("^--data-dir=", "", arg)
    else if (grepl("^--templates-dir=", arg)) config$templates_dir <- sub("^--templates-dir=", "", arg)
    else if (grepl("^--etol=", arg)) config$etol <- sub("^--etol=", "", arg)
    else if (grepl("^--mtol=", arg)) config$mtol <- sub("^--mtol=", "", arg)
    else if (grepl("^--maxe=", arg)) config$maxe <- sub("^--maxe=", "", arg)
    else if (grepl("^--maxm=", arg)) config$maxm <- sub("^--maxm=", "", arg)
    else if (grepl("^--processors=", arg)) config$processors <- sub("^--processors=", "", arg)
    else if (arg %in% c("--help", "-h")) {
      cat(paste(
        "Usage: Rscript WriteFlexmirtFiles.R [options]",
        "  --data-dir=PATH       Generated-data directory",
        "  --templates-dir=PATH  Directory containing the four templates",
        "  --etol=VALUE          E-step tolerance (default: 1e-4)",
        "  --mtol=VALUE          M-step tolerance (default: 1e-7)",
        "  --maxe=N              Maximum E cycles (default: 20000)",
        "  --maxm=N              Maximum M iterations (default: 100)",
        "  --processors=N        Processors per flexMIRT run (default: 1)",
        "  --overwrite           Replace existing generated model files",
        sep = "\n"
      ), "\n")
      quit(save = "no", status = 0L)
    } else stop("Unknown argument: ", arg)
  }
  config$data_dir <- normalizePath(config$data_dir, mustWork = FALSE)
  config$templates_dir <- normalizePath(config$templates_dir, mustWork = FALSE)
  config
}

read_csv <- function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, quote = TRUE, na = "NA")

template_name <- function(qspec, calibration) {
  group_name <- if (calibration == "multigroup") "multigroup" else "singlegroup"
  q_name <- if (qspec == "misspec") "mis" else qspec
  paste0("Flex_template_Q_", q_name, "_", group_name, ".flexmirt")
}

replace_processor_option <- function(lines, processors) {
  matches <- grep("^\\s*Processors\\s*=", lines, ignore.case = TRUE)
  if (length(matches) != 1L) stop("Expected exactly one Processors option")
  lines[matches] <- paste0("  Processors = ", processors, ";")
  lines
}

normalize_template_metadata <- function(lines, qspec, calibration) {
  description <- if (calibration == "multigroup") "multigroup" else "pooled single-group"
  q_label <- if (qspec == "true") "correct Q-matrix" else "misspecified Q-matrix"
  matches <- grep("^\\s*Description\\s*=", lines, ignore.case = TRUE)
  if (length(matches) != 1L) stop("Expected exactly one Description line")
  lines[matches] <- paste0("  Description = \"", description, " calibration with ", q_label, "\";")

  # The source templates contain SavePrm twice with different capitalization.
  save_prm <- grep("^\\s*SavePrm\\s*=", lines, ignore.case = TRUE)
  if (!length(save_prm)) stop("SavePrm option not found")
  if (length(save_prm) > 1L) lines <- lines[-save_prm[-1L]]
  lines
}

grade_block_bounds <- function(lines, grade) {
  start <- grep(paste0("^\\s*%", grade, "%\\s*$"), lines)
  if (length(start) != 1L) stop("Expected one %", grade, "% marker")
  later_markers <- grep("^\\s*%G[345]%\\s*$|^\\s*<Constraints>", lines)
  later_markers <- later_markers[later_markers > start]
  end <- if (length(later_markers)) min(later_markers) - 1L else length(lines)
  c(start = start, end = end)
}

replace_in_grade_block <- function(lines, grade, file_name, n) {
  bounds <- grade_block_bounds(lines, grade)
  indexes <- seq.int(bounds[["start"]], bounds[["end"]])
  file_lines <- indexes[grepl("^\\s*File\\s*=", lines[indexes], ignore.case = TRUE)]
  n_lines <- indexes[grepl("^\\s*N\\s*=", lines[indexes], ignore.case = TRUE)]
  if (length(file_lines) != 1L || length(n_lines) != 1L) {
    stop("Expected one File and one N line in ", grade, " block")
  }
  lines[file_lines] <- paste0("  File = \"", file_name, "\";")
  lines[n_lines] <- paste0("  N = ", n, ";")
  lines
}

read_anchor_items <- function(path) {
  x <- read_csv(path)
  if (!nrow(x)) return(integer(0))
  items <- as.integer(x[[1L]])
  if (anyNA(items) || anyDuplicated(items)) stop("Invalid anchor file: ", path)
  items
}

q_matrix_for <- function(data_dir, qspec) {
  q <- read_csv(file.path(data_dir, paste0("Q_", qspec, ".csv")))
  item_names <- q[[1L]]
  result <- as.matrix(q[, -1L, drop = FALSE])
  storage.mode(result) <- "integer"
  rownames(result) <- item_names
  result
}

effect_constraints <- function(group_left, prefix_left, group_right, prefix_right,
                               anchor_items, q_matrix, label) {
  if (!length(anchor_items)) return(character(0))
  lines <- paste0("  // Generated anchor constraints: ", label)
  for (item in anchor_items) {
    left <- paste0("(", prefix_left, item, ")")
    right <- paste0("(", prefix_right, item, ")")
    lines <- c(lines, sprintf("  Equal %s, %s, Intercept: %s, %s, Intercept;",
                              group_left, left, group_right, right))
    attributes <- which(q_matrix[item, ] == 1L)
    for (attribute in attributes) {
      lines <- c(lines, sprintf("  Equal %s, %s, MainEffect(%d): %s, %s, MainEffect(%d);",
                                group_left, left, attribute, group_right, right, attribute))
    }
    if (length(attributes) >= 2L) {
      pairs <- utils::combn(attributes, 2L)
      for (column in seq_len(ncol(pairs))) {
        a <- pairs[1L, column]; b <- pairs[2L, column]
        lines <- c(lines, sprintf("  Equal %s, %s, Interaction(%d,%d): %s, %s, Interaction(%d,%d);",
                                  group_left, left, a, b, group_right, right, a, b))
      }
    }
  }
  lines
}

insert_after_marker <- function(lines, marker, block) {
  marker_index <- grep(marker, lines, fixed = TRUE)
  if (length(marker_index) != 1L) stop("Expected one marker: ", marker)
  if (!length(block)) return(lines)
  append(lines, block, after = marker_index)
}

create_pooled_responses <- function(rep_dir, n_per_grade, overwrite = FALSE) {
  output_path <- file.path(rep_dir, "responses_all.csv")
  if (file.exists(output_path) && !overwrite) return(output_path)

  grades <- lapply(3:5, function(grade) {
    input <- read_csv(file.path(rep_dir, paste0("responses_G", grade, ".csv")))
    if (nrow(input) != n_per_grade || ncol(input) != 30L) {
      stop("Unexpected response dimensions for Grade ", grade, " in ", rep_dir)
    }
    input
  })
  columns <- c("Grade", paste0("I", 1:30), paste0("J", 1:30), paste0("K", 1:30))
  pooled <- as.data.frame(matrix(-99L, nrow = 3L * n_per_grade, ncol = length(columns)))
  names(pooled) <- columns
  pooled$Grade <- rep(3:5, each = n_per_grade)
  pooled[seq_len(n_per_grade), paste0("I", 1:30)] <- grades[[1L]]
  pooled[n_per_grade + seq_len(n_per_grade), paste0("J", 1:30)] <- grades[[2L]]
  pooled[2L * n_per_grade + seq_len(n_per_grade), paste0("K", 1:30)] <- grades[[3L]]
  write_csv(pooled, output_path)
  output_path
}

validate_model <- function(path, analysis, rep_row, anchors34, anchors45, q_matrix, config) {
  lines <- readLines(path, warn = FALSE)
  validate_tolerances(path, config$etol, config$mtol, config$maxe, config$maxm)
  if (sum(grepl(paste0("^\\s*Processors\\s*=\\s*", config$processors, "\\s*;"), lines, ignore.case = TRUE)) != 1L) {
    stop("Processor validation failed: ", path)
  }
  if (any(grepl("GoogleDrive|\\\\Mac\\\\Home|C001\\\\rep_001", lines))) {
    stop("Stale template path remains in ", path)
  }
  expected_n <- if (analysis$calibration == "multigroup") rep(rep_row$N_per_grade, 3L) else 3L * rep_row$N_per_grade
  actual_n <- as.integer(sub(".*=\\s*([0-9]+)\\s*;.*", "\\1", grep("^\\s*N\\s*=", lines, value = TRUE)))
  if (!identical(actual_n, as.integer(expected_n))) stop("N validation failed: ", path)

  expected_effects <- function(items) sum(vapply(items, function(item) {
    k <- sum(q_matrix[item, ])
    1L + k + choose(k, 2L)
  }, numeric(1)))
  generated <- grep("^\\s*Equal .*", lines, value = TRUE)
  if (length(generated) != expected_effects(anchors34) + expected_effects(anchors45)) {
    stop("Anchor-constraint count failed: ", path)
  }
  TRUE
}

build_one_model <- function(analysis, rep_row, config) {
  rep_dir <- file.path(config$data_dir, rep_row$relative_path)
  if (!dir.exists(rep_dir)) stop("Replication directory not found: ", rep_dir)
  q_matrix <- q_matrix_for(config$data_dir, analysis$Qspec)
  anchors34 <- read_anchor_items(file.path(rep_dir, "anchors_G3_G4.csv"))
  anchors45 <- read_anchor_items(file.path(rep_dir, "anchors_G4_G5.csv"))
  template_path <- file.path(config$templates_dir, template_name(analysis$Qspec, analysis$calibration))
  if (!file.exists(template_path)) stop("Template not found: ", template_path)
  lines <- readLines(template_path, warn = FALSE)
  lines <- normalize_template_metadata(lines, analysis$Qspec, analysis$calibration)
  lines <- replace_processor_option(lines, config$processors)

  if (analysis$calibration == "multigroup") {
    for (grade in 3:5) {
      lines <- replace_in_grade_block(lines, paste0("G", grade),
                                      paste0("responses_G", grade, ".csv"), rep_row$N_per_grade)
    }
    block34 <- effect_constraints("G3", "I", "G4", "I", anchors34, q_matrix, "G3-G4")
    block45 <- effect_constraints("G4", "I", "G5", "I", anchors45, q_matrix, "G4-G5")
    lines <- insert_after_marker(lines, "// Linking attributes between G3 and G4", block34)
    lines <- insert_after_marker(lines, "// Linking attributes between G4 and G5", block45)
  } else {
    create_pooled_responses(rep_dir, rep_row$N_per_grade, config$overwrite)
    lines <- replace_in_grade_block(lines, "G3", "responses_all.csv", 3L * rep_row$N_per_grade)
    block34 <- effect_constraints("G3", "I", "G3", "J", anchors34, q_matrix, "I-J (G3-G4)")
    block45 <- effect_constraints("G3", "J", "G3", "K", anchors45, q_matrix, "J-K (G4-G5)")
    lines <- insert_after_marker(lines, "// Linking attributes between I and J", block34)
    lines <- insert_after_marker(lines, "// Linking attributes between J and K", block45)
  }

  output_name <- sprintf("model_%s_%s_rep%03d.flexmirt",
                         analysis$analysis_condition_id, analysis$data_condition_id, rep_row$rep)
  output_path <- file.path(rep_dir, output_name)
  if (file.exists(output_path) && !config$overwrite) stop("Model already exists: ", output_path)
  writeLines(lines, output_path, useBytes = TRUE)
  update_tolerances(output_path, config$etol, config$mtol, config$maxe, config$maxm)
  validate_model(output_path, analysis, rep_row, anchors34, anchors45, q_matrix, config)

  data.frame(
    analysis_condition_id = analysis$analysis_condition_id,
    data_condition_id = analysis$data_condition_id,
    rep = rep_row$rep,
    Qspec = analysis$Qspec,
    calibration = analysis$calibration,
    model_relative_path = file.path(rep_row$relative_path, output_name),
    response_file = if (analysis$calibration == "multigroup")
      "responses_G3.csv|responses_G4.csv|responses_G5.csv" else "responses_all.csv",
    N_per_grade = rep_row$N_per_grade,
    N_model = if (analysis$calibration == "multigroup") rep_row$N_per_grade else 3L * rep_row$N_per_grade,
    n_anchors_G3_G4 = length(anchors34),
    n_anchors_G4_G5 = length(anchors45),
    MaxE = config$maxe, MaxM = config$maxm, Etol = config$etol, Mtol = config$mtol,
    Processors = config$processors,
    validation_passed = TRUE,
    stringsAsFactors = FALSE
  )
}

main <- function() {
  config <- parse_args(commandArgs(trailingOnly = TRUE))
  analysis_path <- file.path(config$data_dir, "analysis_conditions.csv")
  rep_path <- file.path(config$data_dir, "rep_manifest_all.csv")
  if (!file.exists(analysis_path) || !file.exists(rep_path)) stop("Required manifests not found in ", config$data_dir)
  analyses <- read_csv(analysis_path)
  reps <- read_csv(rep_path)
  if (!all(reps$validation_passed)) stop("Input manifest contains unvalidated replications")

  rows <- vector("list", nrow(reps) * 4L)
  index <- 0L
  for (rep_index in seq_len(nrow(reps))) {
    rep_row <- reps[rep_index, , drop = FALSE]
    selected <- analyses[analyses$data_condition_id == rep_row$data_condition_id, , drop = FALSE]
    if (nrow(selected) != 4L) stop("Expected four analyses for ", rep_row$data_condition_id)
    for (analysis_index in seq_len(nrow(selected))) {
      index <- index + 1L
      rows[[index]] <- build_one_model(selected[analysis_index, , drop = FALSE], rep_row, config)
    }
    if (rep_index %% 25L == 0L || rep_index == nrow(reps)) {
      message(sprintf("Prepared %d/%d datasets (%d model files)", rep_index, nrow(reps), index))
    }
  }
  manifest <- do.call(rbind, rows)
  write_csv(manifest, file.path(config$data_dir, "flexmirt_manifest.csv"))
  message("flexMIRT file generation complete: ", nrow(manifest), " validated model files")
}

if (sys.nframe() == 0L) main()
