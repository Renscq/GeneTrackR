test_that("length extraction distinguishes inclusive spans from spliced lengths", {
  object <- coverage_gene_model()
  before <- serialize(object, NULL)
  genes <- get_gene_length_distribution_table(object, feature = "gene")
  expect_identical(genes$length_bp, c(12, 12, 4))
  spliced <- get_gene_length_distribution_table(object, feature = "transcript")
  genomic <- get_gene_length_distribution_table(object, feature = "transcript",
                                                 transcript_length = "genomic")
  expect_identical(spliced$length_bp, c(9, 9, 4))
  expect_identical(genomic$length_bp, c(12, 12, 4))
  exons <- get_gene_length_distribution_table(object, feature = "exon", unit = "segment")
  expect_identical(exons$length_bp, c(5, 4, 5, 4, 4))
  totals <- get_gene_length_distribution_table(object, feature = "exon", unit = "transcript")
  expect_identical(totals$length_bp, c(9, 9, 4))
  expect_identical(totals$transcript_id, c("plus_tx", "minus_tx", "noncoding_tx"))
  expect_identical(get_genepred_length_table(object, feature = "gene"), genes)
  expect_identical(serialize(object, NULL), before)
})

test_that("CDS and UTR lengths respect exon gaps and strand direction", {
  object <- coverage_gene_model()
  cds <- get_gene_length_distribution_table(object, feature = "cds")
  expect_identical(cds$length_bp, c(5, 5))
  expect_identical(cds$transcript_id, c("plus_tx", "minus_tx"))
  segments <- get_gene_length_distribution_table(object, feature = "cds", unit = "segment")
  expect_identical(segments$start, c(3L, 9L, 23L, 29L))
  expect_identical(segments$end, c(5L, 10L, 25L, 30L))
  expect_identical(segments$length_bp, c(3, 2, 3, 2))
  five <- get_gene_length_distribution_table(object, feature = "five_utr", unit = "segment")
  three <- get_gene_length_distribution_table(object, feature = "three_utr", unit = "segment")
  expect_identical(five$start, c(1L, 31L))
  expect_identical(five$end, c(2L, 32L))
  expect_identical(three$start, c(11L, 21L))
  expect_identical(three$end, c(12L, 22L))
  utr <- get_gene_length_distribution_table(object, feature = "utr", keep_zero = TRUE)
  expect_identical(utr$length_bp, c(4, 4))
  all_features <- get_gene_length_distribution_table(object, feature = "all")
  expect_identical(as.integer(table(all_features$feature)), c(3L, 3L, 5L, 2L, 2L, 2L, 2L))
})

test_that("length queries return empty schemas and reject invalid input", {
  object <- coverage_gene_model()
  for (feature in c("gene", "transcript", "exon", "cds", "utr")) {
    empty <- get_gene_length_distribution_table(object, feature = feature, chrom = "chr9")
    expect_equal(nrow(empty), 0L)
    expect_identical(empty$length_bp, numeric())
    expect_true(all(c("feature", "gene_id", "transcript_id", "start", "end") %in% names(empty)))
  }
  selected <- get_gene_length_distribution_table(object, feature = "gene", chrom = "chr1",
                                                  start = 1L, end = 12L, mode = "within")
  expect_identical(selected$gene_id, "plus_gene")
  expect_error(get_gene_length_distribution_table(Feature(data.frame(
    chrom = "chr1", start = 1L, end = 3L
  ))), "transcript/exon")
  expect_error(get_gene_length_distribution_table(object, feature = "invalid"), "arg")
  expect_error(plot_gene_length_distribution(object, feature = "gene", chrom = "chr9"),
                "No length records")
})

test_that("length plots expose the verified data for every supported plot type", {
  object <- coverage_gene_model()
  for (type in c("density", "histogram", "boxplot", "violin")) {
    observed <- plot_gene_length_distribution(
      object, feature = "gene", group_by = "none", plot_type = type,
      scale = "linear", fill_colors = c(all = "#4477AA"), return_data = TRUE
    )
    expect_s3_class(observed$plot, "ggplot")
    expect_identical(observed$data$length_plot, c(12, 12, 4))
    expect_true(nrow(ggplot2::ggplot_build(observed$plot)$data[[1L]]) > 0L)
  }
  logarithmic <- plot_genepred_length_distribution(object, feature = "all",
    plot_type = "histogram", return_data = TRUE)
  expect_equal(logarithmic$data$length_plot, log10(logarithmic$data$length_bp))
  expect_s3_class(logarithmic$plot$facet, "FacetWrap")
  expect_s3_class(plot_gene_length_distribution(object, feature = "gene"), "ggplot")
})

test_that("multi-exon summaries retain one row per transcript and feature", {
  transcripts <- data.frame(
    transcript_id = c("plus_tx", "minus_tx"), gene_id = c("plus_gene", "minus_gene"),
    chrom = "chr1", strand = c("+", "-"), tx_start = c(1L, 101L), tx_end = c(34L, 134L),
    cds_start = c(13L, 113L), cds_end = c(22L, 122L), exon_count = 4L, gene_type = "coding"
  )
  exons <- data.frame(
    transcript_id = rep(transcripts$transcript_id, each = 4L),
    gene_id = rep(transcripts$gene_id, each = 4L), chrom = "chr1",
    strand = rep(c("+", "-"), each = 4L),
    exon_number = c(1L, 2L, 3L, 4L, 4L, 3L, 2L, 1L),
    exon_start = c(1L, 11L, 21L, 31L, 101L, 111L, 121L, 131L),
    exon_end = c(4L, 14L, 24L, 34L, 104L, 114L, 124L, 134L), exon_frame = NA_integer_
  )
  object <- GenePred(transcripts, exons)
  before <- serialize(object, NULL)
  lengths <- c(exon = 16, cds = 4, utr = 12, five_utr = 6, three_utr = 6)
  starts <- list(exon = c(1L, 101L), cds = c(13L, 113L), utr = c(1L, 101L),
    five_utr = c(1L, 123L), three_utr = c(23L, 101L))
  ends <- list(exon = c(34L, 134L), cds = c(22L, 122L), utr = c(34L, 134L),
    five_utr = c(12L, 134L), three_utr = c(34L, 112L))
  segment_counts <- c(exon = 8L, cds = 4L, utr = 8L, five_utr = 4L, three_utr = 4L)
  for (feature in names(lengths)) {
    observed <- get_gene_length_distribution_table(object, feature = feature, unit = "transcript")
    expect_identical(observed$transcript_id, c("plus_tx", "minus_tx"))
    expect_identical(as.character(observed$feature), rep(feature, 2L))
    expect_identical(observed$unit, rep("transcript", 2L))
    expect_identical(observed$length_bp, rep(unname(lengths[[feature]]), 2L))
    expect_identical(observed$start, starts[[feature]])
    expect_identical(observed$end, ends[[feature]])
    segments <- get_gene_length_distribution_table(object, feature = feature, unit = "segment")
    expect_identical(nrow(segments), unname(segment_counts[[feature]]))
    expect_identical(sum(segments$length_bp), 2 * unname(lengths[[feature]]))
    automatic <- get_gene_length_distribution_table(object, feature = feature)
    expect_identical(automatic, if (feature == "exon") segments else observed)
    expect_identical(get_genepred_length_table(object, feature = feature, unit = "transcript"), observed)
  }
  expect_identical(serialize(object, NULL), before)
})

test_that("transcript aggregation preserves zero filtering, empty schema and input", {
  segments <- get_gene_length_distribution_table(coverage_gene_model(), feature = "cds", unit = "segment")
  segments[, length_bp := c(0, 0, 3, 2)]
  before <- serialize(segments, NULL)
  kept <- GeneTrackR:::summarise_segments_by_transcript(segments, feature = "cds", keep_zero = TRUE)
  expect_identical(kept$transcript_id, c("plus_tx", "minus_tx"))
  expect_identical(kept$length_bp, c(0, 5))
  filtered <- GeneTrackR:::summarise_segments_by_transcript(segments, feature = "cds")
  expect_identical(filtered$transcript_id, "minus_tx")
  expect_identical(filtered$length_bp, 5)
  empty <- GeneTrackR:::summarise_segments_by_transcript(segments[0L], feature = "cds", keep_zero = TRUE)
  expect_identical(names(empty), names(segments))
  expect_identical(empty$length_bp, numeric())
  expect_identical(empty$start, integer())
  expect_identical(serialize(segments, NULL), before)
})

test_that("transcript CDS plots count each transcript once through both APIs", {
  object <- coverage_gene_model()
  for (unit in c("auto", "transcript")) {
    observed <- plot_gene_length_distribution(object, feature = "cds", unit = unit,
      plot_type = "histogram", group_by = "none", scale = "linear", return_data = TRUE)
    expect_identical(observed$data$transcript_id, c("plus_tx", "minus_tx"))
    expect_identical(observed$data$length_plot, c(5, 5))
    expect_equal(sum(ggplot2::ggplot_build(observed$plot)$data[[1L]]$count), 2)
    alias <- plot_genepred_length_distribution(object, feature = "cds", unit = unit,
      plot_type = "histogram", group_by = "none", scale = "linear", return_data = TRUE)
    expect_identical(alias$data, observed$data)
  }
})
