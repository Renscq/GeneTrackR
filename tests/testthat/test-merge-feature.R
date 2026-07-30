make_merge_test_annotation <- function(
    gene_id = "GeneA",
    transcript_id = "TxA1",
    start = 100L,
    source = "test") {
  transcripts <- data.table::data.table(
    transcript_id = transcript_id,
    gene_id = gene_id,
    chrom = "chr1",
    strand = "+",
    tx_start = as.integer(start),
    tx_end = as.integer(start + 400L),
    cds_start = as.integer(start + 50L),
    cds_end = as.integer(start + 350L),
    exon_count = 2L,
    score = NA_real_,
    gene_type = "coding",
    cds_start_stat = "cmpl",
    cds_end_stat = "cmpl",
    exon_frames = "0,0"
  )
  exons <- data.table::data.table(
    transcript_id = rep(transcript_id, 2L),
    gene_id = rep(gene_id, 2L),
    chrom = rep("chr1", 2L),
    strand = rep("+", 2L),
    exon_number = 1:2,
    exon_start = as.integer(c(start, start + 250L)),
    exon_end = as.integer(c(start + 150L, start + 400L)),
    exon_frame = c(0L, 0L)
  )

  GenePred(
    transcripts = transcripts,
    exons = exons,
    meta = list(source_file = source)
  )
}


test_that("merge_feature combines different gene models", {
  gene_a <- make_merge_test_annotation(
    gene_id = "GeneA",
    transcript_id = "TxA1",
    start = 100L
  )
  gene_b <- make_merge_test_annotation(
    gene_id = "GeneB",
    transcript_id = "TxB1",
    start = 1000L
  )

  merged <- merge_feature(
    gene_a,
    gene_b,
    source_names = c("gene_a", "gene_b")
  )

  expect_s3_class(merged, "GenePred")
  expect_setequal(merged$genes$gene_id, c("GeneA", "GeneB"))
  expect_setequal(merged$transcripts$transcript_id, c("TxA1", "TxB1"))
  expect_setequal(merged$data$track_source, c("gene_a", "gene_b"))
})


test_that("merge_feature accepts a named list and mixed Feature inputs", {
  annotation <- make_merge_test_annotation()
  peaks <- Feature(data.table::data.table(
    feature_id = "peak_1",
    name = "peak_1",
    chrom = "chr1",
    start = 700L,
    end = 750L,
    type = "peak",
    strand = "*"
  ))

  merged <- merge_feature(list(annotation = annotation, peaks = peaks))

  expect_s3_class(merged, "Feature")
  expect_true("peak_1" %in% merged$data$feature_id)
  expect_setequal(merged$data$track_source, c("annotation", "peaks"))
})


test_that("merge_feature warns and deduplicates complete transcript models", {
  first <- make_merge_test_annotation(source = "first")
  second <- make_merge_test_annotation(source = "second")

  expect_warning(
    merged <- merge_feature(
      first,
      second,
      source_names = c("first", "second"),
      conflict = "deduplicate"
    ),
    "Duplicated annotation identifiers"
  )

  expect_equal(nrow(merged$genes), 1L)
  expect_equal(nrow(merged$transcripts), 1L)
  expect_equal(nrow(merged$exons), 2L)
  expect_equal(merged$transcripts$track_source, "first")
  expect_false(anyDuplicated(merged$data$feature_id) > 0L)
})


test_that("merge_feature keeps unique transcripts under a duplicated gene", {
  first <- make_merge_test_annotation(
    gene_id = "GeneA",
    transcript_id = "TxA1",
    start = 100L
  )
  second <- make_merge_test_annotation(
    gene_id = "GeneA",
    transcript_id = "TxA2",
    start = 200L
  )

  expect_warning(
    merged <- merge_feature(
      first,
      second,
      source_names = c("first", "second"),
      conflict = "deduplicate"
    ),
    "gene_id=1"
  )

  expect_equal(nrow(merged$genes), 1L)
  expect_setequal(merged$transcripts$transcript_id, c("TxA1", "TxA2"))
  expect_equal(merged$genes$n_transcripts, 2L)
  expect_equal(merged$genes$gene_start, 100L)
  expect_equal(merged$genes$gene_end, 600L)
})


test_that("merge_feature renames conflicts and preserves hierarchy", {
  first <- make_merge_test_annotation(source = "first")
  second <- make_merge_test_annotation(source = "second")

  expect_warning(
    merged <- merge_feature(
      first,
      second,
      source_names = c("first", "second"),
      conflict = "rename"
    ),
    "Applying `conflict = \"rename\"`"
  )

  expect_setequal(
    merged$genes$gene_id,
    c("GeneA", "second_GeneA")
  )
  expect_setequal(
    merged$transcripts$transcript_id,
    c("TxA1", "second_TxA1")
  )
  expect_true(all(
    merged$exons$transcript_id %in%
      merged$transcripts$transcript_id
  ))
  expect_true(all(
    merged$transcripts$gene_id %in%
      merged$genes$gene_id
  ))
  expect_false(anyDuplicated(merged$data$feature_id) > 0L)
})


test_that("merge_feature can reject duplicated identifiers", {
  first <- make_merge_test_annotation()
  second <- make_merge_test_annotation()

  expect_warning(
    expect_error(
      merge_feature(first, second, conflict = "error"),
      "cannot be merged"
    ),
    "Duplicated annotation identifiers"
  )
})
