# Author: Rensc
# Date: 2026-11-21
# Version: dev001
# Function: Regress transcript-specific views of shared annotation features.
# Input: Synthetic GFF3, GRanges and legacy Feature annotations.
# Output: Assertions for membership, provenance and compatibility.

transcript_view_fixture <- function() {
  path <- tempfile(fileext = ".gff3")
  on.exit(unlink(path), add = TRUE)
  writeLines(c(
    "##gff-version 3",
    "chr1\tfixture\tgene\t1\t100\t.\t-\t.\tID=g1",
    "chr1\tfixture\tmRNA\t1\t100\t.\t-\t.\tID=t%2C1;Parent=g1",
    "chr1\tfixture\tmRNA\t2\t100\t.\t-\t.\tID=t2;Parent=g1",
    "chr1\tfixture\tmRNA\t3\t100\t.\t-\t.\tID=t3;Parent=g1",
    "chr1\tfixture\texon\t80\t90\t2.5\t-\t.\tID=shared;Parent=t%2C1,t2,t%2C1",
    "chr1\tfixture\texon\t80\t90\t.\t-\t.\tID=twin;Parent=t%2C1",
    "chr1\tfixture\texon\t50\t60\t.\t-\t.\tID=other;Parent=t2",
    "chr1\tfixture\tCDS\t82\t85\t.\t-\t0\tID=coding;Parent=t2,t%2C1,t2",
    "chr1\tfixture\tCDS\t51\t55\t.\t-\t2\tID=coding;Parent=t2,t%2C1"
  ), path)
  read_gff(path, verbose = FALSE, progress = FALSE)
}

transcript_view <- function(object) as_granges(object, level = "transcript_features")

test_that("shared exons and CDS use one membership per source row and transcript", {
  feature <- transcript_view_fixture()
  original <- lapply(feature, function(value) {
    if (is.data.frame(value)) as.data.frame(value) else value
  })
  flat <- as_granges(feature, level = "feature")
  view <- transcript_view(feature)
  expect_s4_class(view, "GRangesList")
  expect_equal(as_granges(feature, level = c(view = "transcript_features")), view)
  expect_identical(names(view), c("t%2C1", "t2", "t3"))
  expect_identical(unname(lengths(view)), c(4L, 4L, 0L))
  expect_length(flat, 9L)
  mapping <- S4Vectors::mcols(view)$source_rows
  shared <- which(S4Vectors::mcols(flat)$feature_id == "shared")
  expect_identical(sum(mapping[[1L]] == shared), 1L)
  expect_identical(sum(mapping[[2L]] == shared), 1L)
  for (i in seq_along(view)) {
    expect_false(anyDuplicated(mapping[[i]]) > 0L)
    expect_equal(view[[i]], flat[mapping[[i]]])
  }
  expect_identical(as.character(GenomicRanges::strand(view[[1L]])), rep("-", 4L))
  first <- S4Vectors::mcols(view[[1L]])
  expect_identical(first$feature_id, c("coding", "shared", "twin", "coding"))
  expect_identical(first$parent_ids[[2L]], c("t%2C1", "t2", "t%2C1"))
  expect_identical(first$phase[first$type == "CDS"], c("2", "0"))
  expect_equal(lapply(feature, function(value) {
    if (is.data.frame(value)) as.data.frame(value) else value
  }), original)
})

test_that("the view follows feature-level trimming and merged parent identifiers", {
  annotation <- transcript_view_fixture()
  expect_error(retrieve_feature(annotation, level = "feature", chrom = "chr1", start = 81L,
    end = 85L, mode = "trim", as = "Feature"), "discard complete GFF3 parent links")
  flat <- as_granges(annotation, level = "feature")
  region <- as_feature(IRanges::restrict(flat, start = 81L, end = 85L))
  view <- transcript_view(region)
  expect_identical(unname(lengths(view)), c(3L, 2L, 0L))
  expect_true(all(IRanges::start(unlist(view)) >= 81L & IRanges::end(unlist(view)) <= 85L))
  expect_equal(view[[1L]], as_granges(region, level = "feature")[S4Vectors::mcols(view)$source_rows[[1L]]])

  expect_warning(merged <- merge_feature(annotation, annotation,
    source_names = c("A", "B"), conflict = "rename"), "Duplicated annotation")
  merged_view <- transcript_view(merged)
  expect_identical(unname(lengths(merged_view[c("B_t%2C1", "B_t2")])), c(4L, 4L))
  expect_identical(sum(S4Vectors::mcols(merged_view[["B_t2"]])$feature_id == "B_shared"), 1L)
})

test_that("legacy single-parent annotations and default projections stay compatible", {
  annotation <- transcript_view_fixture()
  legacy <- annotation
  legacy$data <- data.table::copy(annotation$data)
  legacy$data$parent_ids <- NULL
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(legacy, path)
  restored <- readRDS(path)
  view <- transcript_view(restored)
  expect_identical(unname(lengths(view)), c(2L, 3L, 0L))
  expect_false("shared" %in% S4Vectors::mcols(view[["t2"]])$feature_id)
  expect_equal(as_granges(annotation), as_granges(annotation, level = "gene"))
  expect_equal(as_granges(annotation, level = "trans"), as_granges(restored, level = "transcript"))
  for (level in c("gene", "transcript", "exon")) {
    expect_equal(summary_feature(annotation, level = level), summary_feature(restored, level = level))
  }

  flat <- FeatureTrack(data.frame(chrom = "chr1", start = c(3L, 5L), end = c(8L, 6L),
    type = c("feature", "CDS"), level = c("exon", "subfeature"), strand = "+",
    transcript_id = c("old_tx", NA_character_), parent_id = c(NA_character_, "old_tx")))
  class(flat) <- "FeatureTrack"
  expect_identical(names(transcript_view(flat)), "old_tx")
  expect_identical(length(transcript_view(flat)[[1L]]), 2L)

  legacy$data$level <- NULL
  expect_identical(unname(lengths(transcript_view(legacy))), c(2L, 3L, 0L))
})

test_that("single-parent fields use their respective identifier namespaces", {
  feature <- Feature(data.frame(chrom = "chr1", start = c(1L, 1L, 4L, 6L),
    end = c(20L, 20L, 8L, 9L), strand = "+", type = c("mRNA", "mRNA", "exon", "exon"),
    feature_id = c("raw", "one", "parent_only", "explicit_transcript"),
    transcript_id = c("one", "two", NA_character_, "one"),
    parent_id = c(NA_character_, NA_character_, "one", "one")))
  grouped <- transcript_view(feature)
  expect_identical(S4Vectors::mcols(grouped[["one"]])$feature_id, "explicit_transcript")
  expect_identical(S4Vectors::mcols(grouped[["two"]])$feature_id, "parent_only")

  feature$meta$format <- "GFF"
  grouped <- transcript_view(feature)
  expect_length(grouped[["one"]], 0L)
  expect_identical(S4Vectors::mcols(grouped[["two"]])$feature_id,
    c("parent_only", "explicit_transcript"))
})

test_that("common GFF3 RNA types can declare members of the grouped view", {
  rows <- data.frame(chrom = "chr1", start = c(1L, 2L, 3L, 5L), end = c(20L, 20L, 20L, 10L),
    type = c("snoRNA", "snRNA", "pseudogenic_transcript", "exon"),
    feature_id = c("sno", "sn", "pseudo", "shared"), strand = "+")
  rows$parent_ids <- list(character(), character(), character(), c("sno", "sn", "pseudo"))
  view <- transcript_view(Feature(rows))
  expect_identical(names(view), c("sno", "sn", "pseudo"))
  expect_identical(unname(lengths(view)), rep(1L, 3L))
})

test_that("discontinuous declarations for one transcript form a single group", {
  path <- tempfile(fileext = ".gff3")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("##gff-version 3",
    "chr1\tfixture\tgene\t1\t80\t.\t+\t.\tID=g1",
    "chr1\tfixture\tmRNA\t1\t20\t.\t+\t.\tID=t1;Parent=g1",
    "chr1\tfixture\tmRNA\t60\t80\t.\t+\t.\tID=t1;Parent=g1",
    "chr1\tfixture\texon\t1\t20\t.\t+\t.\tID=e1;Parent=t1",
    "chr1\tfixture\texon\t60\t80\t.\t+\t.\tID=e2;Parent=t1"), path)
  annotation <- read_gff(path, verbose = FALSE, progress = FALSE)
  view <- transcript_view(annotation)
  expect_identical(names(view), "t1")
  expect_identical(IRanges::start(view[["t1"]]), c(1L, 60L))
  expect_length(as_granges(annotation, level = "feature"), 5L)
})

test_that("GFF3 SO accessions and percent-escaped identifiers retain source metadata", {
  path <- tempfile(fileext = ".gff3")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("##gff-version 3",
    "chr1\tfixture\tSO:0000704\t1\t20\t.\t+\t.\tID=g1",
    "chr1\tfixture\tSO:0000234\t1\t20\t.\t+\t.\tID=t%2C1;Parent=g1",
    "chr1\tfixture\tSO:0000147\t1\t8\t.\t+\t.\tID=e1;Parent=t%2c1,t%2C1",
    "chr1\tfixture\tSO:0000316\t3\t8\t.\t+\t0\tID=c1;Parent=t%2c1"), path)
  annotation <- read_gff(path, verbose = FALSE, progress = FALSE)
  view <- transcript_view(annotation)
  expect_identical(names(view), "t%2C1")
  expect_identical(length(view[[1L]]), 2L)
  metadata <- S4Vectors::mcols(view[[1L]])
  expect_identical(metadata$type, c("SO:0000147", "SO:0000316"))
  expect_identical(metadata$parent_ids[[1L]], c("t%2c1", "t%2C1"))
})

test_that("standard ranges retain their identity and metadata inside groups", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1L, 5L), c(20L, 8L)),
    strand = "+", seqlengths = c(chr1 = 100L),
    feature_id = c("raw_tx", "exon1"), type = c("mRNA", "exon"),
    transcript_id = c("tx_id", "tx_id"), parent_ids = IRanges::CharacterList(character(), "raw_tx"),
    source_rows = c("existing1", "existing2"))
  names(ranges) <- c("transcript range", "exon range")
  S4Vectors::metadata(ranges)$genome_build <- "fixture"
  feature <- as_feature(ranges)
  grouped <- transcript_view(feature)
  expect_identical(names(grouped), "tx_id")
  expect_identical(names(grouped[[1L]]), "exon range")
  expect_identical(S4Vectors::mcols(grouped[[1L]])$source_rows, "existing2")
  expect_identical(S4Vectors::mcols(grouped)$source_rows[[1L]], 2L)
  expect_identical(GenomicRanges::seqinfo(grouped), GenomicRanges::seqinfo(ranges))
  expect_identical(S4Vectors::metadata(grouped)$genome_build, "fixture")
})

test_that("unresolvable and ambiguous transcript annotations fail explicitly", {
  annotation <- transcript_view_fixture()
  child <- which(annotation$data$feature_id == "shared")
  edited <- annotation
  edited$data <- data.table::copy(annotation$data)
  edited$data$parent_ids[[child]] <- c("t%2C1", "unknown_tx")
  expect_error(transcript_view(edited), "Unknown.*unknown_tx")
  edited$data$parent_ids[[child]] <- c("t%2C1", NA_character_)
  expect_error(transcript_view(edited), "parent_ids.*nonempty")
  edited$data$parent_ids[child] <- list(character())
  edited$data$parent_id[child] <- NA_character_
  edited$data$transcript_id[child] <- NA_character_
  expect_error(transcript_view(edited), "Missing transcript.*shared")

  tx_row <- which(annotation$data$feature_id == "t2")
  edited$data <- data.table::rbindlist(list(annotation$data, annotation$data[tx_row]))
  edited$data$transcript_id[nrow(edited$data)] <- "conflicting_id"
  expect_error(transcript_view(edited), "Ambiguous transcript.*t2")
  edited$data <- data.table::copy(annotation$data)
  edited$data$strand[child] <- "+"
  expect_error(transcript_view(edited), "incompatible chromosome or strand")
  unlabeled <- Feature(data.frame(chrom = "chr1", start = c(1L, 5L), end = c(10L, 8L),
    type = c("mRNA", "exon"), feature_id = c(NA_character_, "unlinked"),
    transcript_id = c("tx", NA_character_)))
  expect_error(transcript_view(unlabeled), "Missing transcript.*unlinked")
  expect_error(transcript_view(data.frame(chrom = "chr1", pos = 10L)), "must be a GenePred")
  variant <- VariantTrack(data.frame(chrom = "chr1", pos = 10L))
  expect_equal(transcript_view(variant), as_granges(variant))
})

test_that("transcript groups handle empty annotations and older GenePred data", {
  empty <- Feature(data.frame(chrom = character(), start = integer(), end = integer()))
  empty_view <- transcript_view(empty)
  expect_s4_class(empty_view, "GRangesList")
  expect_length(empty_view, 0L)

  path <- system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR", mustWork = TRUE)
  annotation <- read_genepred(path, format = "genePredExt", verbose = FALSE, progress = FALSE)
  expected <- transcript_view(annotation)
  annotation$data <- NULL
  expect_equal(transcript_view(annotation), expected)
})
