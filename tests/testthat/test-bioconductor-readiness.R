test_that("Bioconductor metadata baseline is present", {
  desc <- utils::packageDescription("GeneTrackR")

  expect_identical(desc[["NeedsCompilation"]], "no")
  expect_true(is.null(desc[["Remotes"]]))
  expect_true(is.null(desc[["LinkingTo"]]))

  imports <- trimws(strsplit(desc[["Imports"]], ",", fixed = TRUE)[[1L]])
  expect_false("Rcpp" %in% imports)

  views <- trimws(strsplit(desc[["biocViews"]], ",", fixed = TRUE)[[1L]])
  expect_true(all(c(
    "Software",
    "GenomeAnnotation",
    "Coverage",
    "Visualization"
  ) %in% views))
  expect_true(length(setdiff(views, "Software")) >= 2L)

  description <- desc[["Description"]]
  sentences <- strsplit(description, "[.!?]+")[[1L]]
  sentences <- sentences[nzchar(trimws(sentences))]
  expect_true(length(sentences) >= 3L)
})

test_that("annotation objects interoperate with GenomicRanges", {
  gp <- read_genepred(
    gtr_extdata("gtr_demo.genePredExt"),
    format = "genePredExt",
    verbose = FALSE,
    progress = FALSE
  )
  feature <- read_bed(
    gtr_extdata("gtr_demo_features.bed"),
    verbose = FALSE,
    progress = FALSE
  )

  expect_s4_class(as_granges(gp, level = "gene"), "GRanges")
  expect_s4_class(as_granges(gp, level = "transcript"), "GRanges")
  expect_s4_class(as_granges(feature, level = "feature"), "GRanges")
})
