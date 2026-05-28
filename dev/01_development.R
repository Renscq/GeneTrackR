# Author: Rensc
# Date: 2026-05-26
# Version: 0.1.0
# Function: Development workflow for GeneTrackR
# Input: R package source directory
# Output: Updated documentation, namespace, checks, and build artifacts

# Install development dependencies if needed.
# install.packages(c("devtools", "roxygen2", "testthat"))

library(devtools)
library(roxygen2)

# Document package and update NAMESPACE/man files.
devtools::document()

# Run package checks during active development.
devtools::check(document = FALSE)

# Install package locally.
devtools::install(build_vignettes = FALSE, upgrade = "never")
