test_that("GRanges records keep metadata and names aligned through Feature sorting", {
  sequence_info <- methods::new(
    "Seqinfo",
    seqnames = c("chr1", "chr2", "chrUnused"),
    seqlengths = c(100L, 200L, 300L),
    is_circular = c(FALSE, TRUE, NA),
    genome = rep("synthetic", 3L)
  )
  ranges <- GenomicRanges::GRanges(
    seqnames = c("chr2", "chr1", "chr1"),
    ranges = IRanges::IRanges(c(30L, 7L, 7L), c(35L, 9L, 9L)),
    strand = c("-", "+", "*"),
    seqinfo = sequence_info
  )
  names(ranges) <- c("range_c", "range_a", "range_b")
  metadata <- S4Vectors::make_zero_col_DFrame(3L)
  metadata$feature_id <- c("third", "duplicate", "duplicate")
  metadata$name <- c("label_c", "label_a", "label_b")
  metadata$score <- c(3.5, NA_real_, 1.5)
  metadata$category <- factor(c("c", "a", "b"), levels = c("a", "b", "c", "unused"))
  metadata$members <- IRanges::IntegerList(3L, c(1L, 2L), integer())
  metadata$quality <- S4Vectors::Rle(c(30L, 10L, 20L))
  metadata$details <- list(list(label = "c"), list(label = "a"), NULL)
  S4Vectors::mcols(ranges) <- metadata
  S4Vectors::metadata(ranges) <- list(source = "synthetic", parameters = list(threshold = 2))
  original <- serialize(ranges, NULL)

  feature <- as_feature(ranges)
  observed <- as_granges(feature, level = "feature")

  expect_s3_class(feature, "Feature")
  expect_s3_class(feature, "FeatureTrack")
  expect_identical(feature$data$chrom, c("chr1", "chr1", "chr2"))
  expect_identical(feature$data$feature_id, c("duplicate", "duplicate", "third"))
  expect_identical(as.integer(IRanges::start(observed)), c(7L, 7L, 30L))
  expect_identical(as.character(GenomicRanges::strand(observed)), c("+", "*", "-"))
  expect_identical(names(observed), c("range_a", "range_b", "range_c"))
  expect_identical(S4Vectors::mcols(observed)$name, c("label_a", "label_b", "label_c"))
  expect_identical(S4Vectors::mcols(observed)$score, c(NA_real_, 1.5, 3.5))
  expect_identical(S4Vectors::mcols(observed)$category,
                   factor(c("a", "b", "c"), levels = c("a", "b", "c", "unused")))
  expect_identical(S4Vectors::mcols(observed)$members,
                   IRanges::IntegerList(c(1L, 2L), integer(), 3L))
  expect_identical(S4Vectors::mcols(observed)$quality, S4Vectors::Rle(c(10L, 20L, 30L)))
  expect_identical(S4Vectors::mcols(observed)$details,
                   list(list(label = "a"), NULL, list(label = "c")))
  expect_identical(GenomicRanges::seqinfo(observed), sequence_info)
  expect_identical(S4Vectors::metadata(observed),
                   list(source = "synthetic", parameters = list(threshold = 2)))
  expect_identical(serialize(ranges, NULL), original)
})

test_that("feature-level projection keeps all native Feature metadata for a round trip", {
  feature <- Feature(data.frame(
    chrom = c("chr1", "chr1"), start = c(20L, 5L), end = c(25L, 9L),
    feature_id = c("b", "a"), name = c("B", "A"), type = "peak",
    score = c(4.5, 2.5), source = c("assay_b", "assay_a"),
    parent_id = c("parent_b", "parent_a"), phase = c("2", "0"),
    confidence = c(0.8, 0.9)
  ))
  ranges <- as_granges(feature, level = "feature")
  restored <- as_feature(ranges)

  expect_identical(S4Vectors::mcols(ranges)$score, c(2.5, 4.5))
  expect_identical(S4Vectors::mcols(ranges)$source, c("assay_a", "assay_b"))
  expect_identical(S4Vectors::mcols(ranges)$phase, c("0", "2"))
  expect_identical(S4Vectors::mcols(ranges)$confidence, c(0.9, 0.8))
  expect_true(all(names(feature$data) %in% names(restored$data)))
  for (column in names(feature$data)) {
    expect_identical(restored$data[[column]], feature$data[[column]], info = column)
  }
})

test_that("Feature subsets retain GRanges context without stale row metadata", {
  ranges <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(c(20L, 5L), c(25L, 9L)),
    seqlengths = c(chr1 = 100L)
  )
  names(ranges) <- c("late", "early")
  S4Vectors::mcols(ranges)$category <- c("late_group", "early_group")
  S4Vectors::metadata(ranges) <- list(source = "subset_fixture")
  feature <- as_feature(ranges)
  selected <- retrieve_feature(feature, chrom = "chr1", start = 1L, end = 10L)
  selected$data$category <- "edited_group"
  observed <- as_granges(selected, level = "feature")

  expect_length(observed, 1L)
  expect_identical(names(observed), "early")
  expect_identical(S4Vectors::mcols(observed)$category, "edited_group")
  expect_identical(GenomicRanges::seqinfo(observed), GenomicRanges::seqinfo(ranges))
  expect_identical(S4Vectors::metadata(observed), list(source = "subset_fixture"))
  expect_identical(as_granges(feature, level = "feature")$category,
                   c("early_group", "late_group"))

  empty <- retrieve_feature(feature, chrom = "chr1", start = 50L, end = 60L)
  empty_ranges <- as_granges(empty, level = "feature")
  expect_length(empty_ranges, 0L)
  expect_identical(GenomicRanges::seqinfo(empty_ranges), GenomicRanges::seqinfo(ranges))
  expect_identical(S4Vectors::metadata(empty_ranges), list(source = "subset_fixture"))
})

test_that("empty and unnamed GRanges remain valid feature-level inputs", {
  ranges <- GenomicRanges::GRanges()
  S4Vectors::mcols(ranges)$category <- factor(character(), levels = c("a", "b"))
  feature <- as_feature(ranges)
  observed <- as_granges(feature, level = "feature")
  expect_equal(nrow(feature$data), 0L)
  expect_length(observed, 0L)
  expect_null(names(observed))
  expect_identical(S4Vectors::mcols(observed)$category,
                   factor(character(), levels = c("a", "b")))
})

test_that("ordinary GRanges enter plotting and overlap workflows without invented gene models", {
  ranges <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(c(20L, 5L), c(25L, 9L))
  )
  names(ranges) <- c("late", "early")
  feature <- as_feature(ranges)
  figure <- plot_feature_track(feature, chrom = "chr1", start = 1L, end = 30L)
  expect_s3_class(figure, "ggplot")
  expect_identical(sort(figure$data$start), c(5L, 20L))
  expect_equal(nrow(feature$transcripts), 0L)
  expect_equal(nrow(feature$exons), 0L)
  expect_error(as_genepred(feature), "transcript-level")

  observed <- as_granges(feature, level = "feature")
  query <- GenomicRanges::GRanges("chr1", IRanges::IRanges(8L, 8L))
  hits <- GenomicRanges::findOverlaps(observed, query)
  expect_identical(as.integer(S4Vectors::queryHits(hits)), 1L)
})

test_that("GRanges conversion rejects unrepresentable coordinates and ambiguous columns", {
  zero_width <- GenomicRanges::GRanges("chr1", IRanges::IRanges(5L, width = 0L))
  zero_start <- GenomicRanges::GRanges("chr1", IRanges::IRanges(0L, 3L))
  expect_error(as_feature(zero_width), "positive-width|positive width")
  expect_error(as_feature(zero_start), "1-based|positive")

  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 3L))
  S4Vectors::mcols(ranges)$chrom <- "metadata_chromosome"
  expect_error(as_feature(ranges), "chrom|structural")
  S4Vectors::mcols(ranges) <- NULL
  S4Vectors::mcols(ranges)$score <- "not_numeric"
  expect_error(as_feature(ranges), "score.*numeric")
  S4Vectors::mcols(ranges) <- NULL
  S4Vectors::mcols(ranges)$matrix_value <- matrix(c(1L, 2L), nrow = 1L)
  expect_error(as_feature(ranges), "matrix_value|one-dimensional")
  S4Vectors::mcols(ranges) <- NULL
  S4Vectors::mcols(ranges)$exon_number <- 2.5
  expect_error(as_feature(ranges), "exon_number.*integers")
})

test_that("source-row bookkeeping cannot overwrite an input metadata column", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(20L, 5L), c(25L, 9L)))
  names(ranges) <- c("late", "early")
  S4Vectors::mcols(ranges)[[".gtr_granges_names"]] <- c("metadata_late", "metadata_early")
  S4Vectors::mcols(ranges)[[".gtr_granges_row_id"]] <- c(200L, 100L)
  observed <- as_granges(as_feature(ranges), level = "feature")
  expect_identical(names(observed), c("early", "late"))
  expect_identical(S4Vectors::mcols(observed)[[".gtr_granges_names"]],
                   c("metadata_early", "metadata_late"))
  expect_identical(S4Vectors::mcols(observed)[[".gtr_granges_row_id"]], c(100L, 200L))
})

test_that("missing numeric fields normalize safely and overflowing exon indices fail without warnings", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 3L))
  S4Vectors::mcols(ranges)$score <- NA_character_
  S4Vectors::mcols(ranges)$exon_number <- NA_character_
  expect_no_warning(feature <- as_feature(ranges))
  expect_identical(feature$data$score, NA_real_)
  expect_identical(feature$data$exon_number, NA_integer_)
  S4Vectors::mcols(ranges)$exon_number <- 2^40
  expect_no_warning(expect_error(as_feature(ranges), "exon_number.*integers"))
})

test_that("explicit exon and CDS records support gene models without altering feature metadata", {
  ranges <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 4L),
    ranges = IRanges::IRanges(c(1L, 1L, 1L, 3L), c(12L, 12L, 12L, 8L)),
    strand = "+"
  )
  metadata <- S4Vectors::make_zero_col_DFrame(4L)
  metadata$feature_id <- c("gene", "transcript", "exon", "cds")
  metadata$type <- c("gene", "transcript", "exon", "CDS")
  metadata$gene_id <- rep("gene", 4L)
  metadata$transcript_id <- c(NA_character_, "transcript", "transcript", "transcript")
  metadata$phase <- c(NA_character_, NA_character_, NA_character_, "0")
  metadata$confidence <- c(0.1, 0.2, 0.3, 0.4)
  S4Vectors::mcols(ranges) <- metadata
  feature <- as_feature(ranges)
  model <- as_genepred(feature)
  expect_identical(model$transcripts$cds_start, 3L)
  expect_identical(model$transcripts$cds_end, 8L)
  expect_identical(model$exons$exon_start, 1L)
  expect_identical(model$exons$exon_end, 12L)
  observed <- as_granges(feature, level = "feature")
  observed_metadata <- S4Vectors::mcols(observed)
  expect_identical(observed_metadata$confidence[match(c("gene", "transcript", "exon", "cds"), observed_metadata$feature_id)],
                   c(0.1, 0.2, 0.3, 0.4))

  rebuilt <- retrieve_feature(feature, gene_id = "gene")
  expect_error(as_granges(rebuilt, level = "feature"), "metadata were lost")
})

test_that("numeric annotation text follows Feature normalization without losing valid values", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(20L, 5L), c(25L, 9L)))
  S4Vectors::mcols(ranges)$score <- c("4.5", "2.5")
  S4Vectors::mcols(ranges)$exon_number <- c("2", "1")
  feature <- as_feature(ranges)
  expect_identical(feature$data$score, c(2.5, 4.5))
  expect_identical(feature$data$exon_number, c(1L, 2L))
})

test_that("native Feature export rejects duplicate columns instead of dropping metadata", {
  data <- data.table::data.table(
    chrom = "chr1", start = 1L, end = 3L,
    duplicate = 1L, duplicate = 2L
  )
  feature <- Feature(data)
  expect_error(as_granges(feature, level = "feature"), "unique|duplicat")
})

test_that("native Feature export never returns invalid reserved GRanges metadata", {
  for (column in c("width", "seqnames")) {
    data <- data.frame(chrom = "chr1", start = 1L, end = 3L)
    data[[column]] <- 1L
    expect_error(as_granges(Feature(data), level = "feature"), column)
  }
})

test_that("annotated metadata columns retain names and nested annotations after sorting and subsets", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(20L, 5L), c(25L, 9L)))
  members <- IRanges::IntegerList(red = 2L, blue = 1L)
  S4Vectors::mcols(members)$code <- c("late_code", "early_code")
  S4Vectors::metadata(members) <- list(method = "observed")
  columns <- S4Vectors::make_zero_col_DFrame(2L)
  columns$members <- members
  columns$measurement <- stats::setNames(c(4.5, 2.5), c("late_name", "early_name"))
  S4Vectors::mcols(columns)$description <- c("membership", "measurement")
  S4Vectors::metadata(columns) <- list(dictionary = "fixture")
  S4Vectors::mcols(ranges) <- columns
  feature <- as_feature(ranges)
  observed <- as_granges(feature, level = "feature")
  observed_columns <- S4Vectors::mcols(observed)
  expect_identical(names(observed_columns$members), c("blue", "red"))
  expect_identical(S4Vectors::mcols(observed_columns$members)$code,
                   c("early_code", "late_code"))
  expect_identical(S4Vectors::metadata(observed_columns$members), list(method = "observed"))
  expect_identical(names(observed_columns$measurement), c("early_name", "late_name"))
  expect_identical(S4Vectors::metadata(observed_columns), list(dictionary = "fixture"))
  annotations <- S4Vectors::mcols(observed_columns)
  expect_identical(annotations$description[match(c("members", "measurement"), names(observed_columns))],
                   c("membership", "measurement"))
  selected <- retrieve_feature(feature, chrom = "chr1", start = 1L, end = 10L)
  selected$data$members[[1L]] <- c(7L, 8L)
  selected_ranges <- as_granges(selected, level = "feature")
  selected_members <- S4Vectors::mcols(selected_ranges)$members
  expect_identical(as.list(selected_members), list(blue = c(7L, 8L)))
  expect_identical(S4Vectors::mcols(selected_members)$code, "early_code")
  expect_identical(S4Vectors::metadata(selected_members), list(method = "observed"))
})

test_that("unnamed GRanges provenance detects lossy gene-model reconstruction", {
  ranges <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(1L, 1L, 1L, 3L), c(12L, 12L, 12L, 8L)),
    strand = "+"
  )
  S4Vectors::mcols(ranges) <- S4Vectors::DataFrame(
    feature_id = c("gene", "transcript", "exon", "cds"),
    type = c("gene", "transcript", "exon", "CDS"),
    gene_id = rep("gene", 4L),
    transcript_id = c(NA_character_, "transcript", "transcript", "transcript"),
    score = as.numeric(1:4)
  )
  feature <- as_feature(ranges)
  rebuilt <- retrieve_feature(feature, gene_id = "gene")
  expect_error(as_granges(rebuilt, level = "feature"), "metadata were lost")
})

test_that("custom column names cannot mask Feature construction or interval filters", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(20L, 5L), c(25L, 9L)))
  names(ranges) <- c("late", "early")
  S4Vectors::mcols(ranges)$gene_id <- c("late_gene", "early_gene")
  for (column in c("x", "dt", "keep", "hit", "values", "field", "query_start", "query_end", "start_col", "end_col")) {
    S4Vectors::mcols(ranges)[[column]] <- c(999L, 888L)
  }
  feature <- as_feature(ranges)
  selected <- retrieve_feature(feature, chrom = "chr1", start = 1L, end = 10L,
                                gene_id = "early_gene", pattern = "early")
  observed <- as_granges(selected, level = "feature")
  expect_identical(as.integer(IRanges::start(observed)), 5L)
  expect_identical(names(observed), "early")
  expect_identical(S4Vectors::mcols(observed)$x, 888L)
  expect_identical(S4Vectors::mcols(observed)$dt, 888L)
})
