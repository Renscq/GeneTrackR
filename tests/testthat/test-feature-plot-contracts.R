feature_plot_fixture <- function() {
  Feature(data.frame(
    chrom = rep("chr1", 3L),
    start = c(31L, 11L, 21L),
    end = c(35L, 20L, 23L),
    name = c("outside", "overlap", "contained"),
    strand = c("+", "+", "-")
  ))
}

test_that("feature rectangles keep source bounds and center their labels", {
  track <- feature_plot_fixture()
  strand_colors <- setNames(c("#245A81", "#B65D0D"), c("+", "-"))

  overlapped <- plot_feature_track(
    track,
    chrom = "chr1",
    start = 15L,
    end = 25L,
    mode = "overlap",
    color_by = "strand",
    feature_colors = strand_colors,
    label_by = "name"
  )
  layers <- ggplot2::ggplot_build(overlapped)$data
  rectangles <- layers[[1L]]
  labels <- layers[[2L]]

  expect_equal(overlapped$coordinates$limits$x, c(15L, 25L))
  expect_equal(rectangles$xmin, c(11L, 21L))
  expect_equal(rectangles$xmax, c(20L, 23L))
  expect_identical(rectangles$fill, unname(strand_colors))
  expect_identical(labels$label, c("overlap", "contained"))
  expect_equal(labels$x, c(15.5, 22))
  expect_equal(
    as.numeric(labels$y),
    (rectangles$ymin + rectangles$ymax) / 2
  )

  contained <- plot_feature_track(
    track,
    chrom = "chr1",
    start = 15L,
    end = 25L,
    mode = "within"
  )
  expect_identical(contained$data$name, "contained")
  expect_equal(contained$data$start, 21L)
  expect_equal(contained$data$end, 23L)
})

test_that("empty feature panels preserve the requested interval and draw no rectangles", {
  track <- feature_plot_fixture()
  empty <- plot_feature_track(
    track,
    chrom = "chr2",
    start = 100L,
    end = 120L,
    label_by = "name"
  )

  expect_equal(empty$coordinates$limits$x, c(100L, 120L))
  expect_identical(empty$labels$x, "Chromosome chr2 position (bp)")
  expect_length(empty$layers, 1L)
  expect_s3_class(empty$layers[[1L]]$geom, "GeomBlank")
  expect_error(
    plot_feature_track(track, "chr1", 15L, 25L, color_by = "absent_group"),
    "column was not found in feature data: absent_group"
  )
})
