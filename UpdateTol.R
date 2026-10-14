## Update and validate flexMIRT numerical settings.
## This file has no package dependencies and can be sourced by other scripts.

replace_flex_option <- function(lines, option, value) {
  pattern <- paste0("^\\s*", option, "\\s*=")
  matches <- grep(pattern, lines, ignore.case = TRUE)
  if (length(matches) != 1L) {
    stop("Expected exactly one ", option, " option; found ", length(matches))
  }
  indent <- sub("^(\\s*).*", "\\1", lines[matches])
  lines[matches] <- paste0(indent, option, " = ", value, ";")
  lines
}

update_tolerances <- function(file_path, etol = "1e-4", mtol = "1e-7",
                              maxe = "20000", maxm = "100") {
  if (!file.exists(file_path)) stop("flexMIRT file not found: ", file_path)
  lines <- readLines(file_path, warn = FALSE)
  lines <- replace_flex_option(lines, "MaxE", maxe)
  lines <- replace_flex_option(lines, "MaxM", maxm)
  lines <- replace_flex_option(lines, "Etol", etol)
  lines <- replace_flex_option(lines, "Mtol", mtol)
  writeLines(lines, file_path, useBytes = TRUE)
  invisible(validate_tolerances(file_path, etol, mtol, maxe, maxm))
}

validate_tolerances <- function(file_path, etol = "1e-4", mtol = "1e-7",
                                maxe = "20000", maxm = "100") {
  lines <- readLines(file_path, warn = FALSE)
  expected <- c(MaxE = maxe, MaxM = maxm, Etol = etol, Mtol = mtol)
  for (option in names(expected)) {
    pattern <- paste0("^\\s*", option, "\\s*=\\s*", expected[[option]], "\\s*;")
    if (sum(grepl(pattern, lines, ignore.case = TRUE)) != 1L) {
      stop("Tolerance validation failed for ", option, " in ", file_path)
    }
  }
  TRUE
}

