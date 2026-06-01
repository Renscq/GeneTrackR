gtr_extdata <- function(filename) {
  system.file("extdata", filename, package = "GeneTrackR", mustWork = TRUE)
}

expect_gt_plot <- function(x) {
  testthat::expect_true(
    inherits(x, "ggplot") || inherits(x, "patchwork") || inherits(x, "GeneTrackRPhenoPlot")
  )
}
