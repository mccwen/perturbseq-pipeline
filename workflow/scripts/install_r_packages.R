#!/usr/bin/env Rscript

package <- "mixtools"
required_version <- "2.0.0.1"
installed_version <- if (requireNamespace(package, quietly = TRUE)) {
  as.character(packageVersion(package))
} else {
  NA_character_
}

if (is.na(installed_version) || installed_version != required_version) {
  remotes::install_version(
    package,
    version = required_version,
    repos = "https://cloud.r-project.org",
    dependencies = FALSE,
    upgrade = "never",
    quiet = TRUE
  )
}

if (!requireNamespace(package, quietly = TRUE) ||
    as.character(packageVersion(package)) != required_version) {
  stop("Required R package version is unavailable: mixtools ", required_version)
}

message("Required CRAN package is installed: mixtools ", required_version)
