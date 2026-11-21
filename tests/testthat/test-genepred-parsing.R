# Author: Rensc
# Date: 2026-11-21
# Version: dev001
# Function: Regress public GenePred parsing and validation.
# Input: Synthetic GenePred and GenePredExt records.
# Output: Assertions for coordinates, identifiers, exon frames and invalid records.

test_that("standard GenePred detects the format and converts UCSC intervals", {
  path <- tempfile(fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(c(
    "plus_tx\tchr1\t+\t0\t40\t30\t39\t2\t0,30,\t1,40,",
    "minus_tx\tchr2\t-\t50\t92\t92\t92\t3\t50,70,86\t59,78,92"
  ), path)

  annotation <- read_genepred(path, verbose = FALSE, progress = FALSE)
  transcripts <- as_transcript_table(annotation)
  exons <- as_exon_table(annotation)

  expect_identical(class(annotation), c("Feature", "FeatureTrack", "GenePred"))
  expect_identical(annotation$meta$format, "genePred")
  expect_identical(annotation$meta$coordinate_input, "ucsc")
  expect_identical(annotation$meta$coordinate_internal, "1-based closed")
  expect_identical(annotation$meta$gene_col, "name")
  expect_identical(transcripts$transcript_id, c("plus_tx", "minus_tx"))
  expect_identical(transcripts$gene_id, transcripts$transcript_id)
  expect_identical(transcripts$tx_start, c(1L, 51L))
  expect_identical(transcripts$tx_end, c(40L, 92L))
  expect_identical(transcripts$cds_start[1L], 31L)
  expect_identical(transcripts$cds_end[1L], 39L)
  expect_identical(transcripts$gene_type, c("coding", "non-coding"))
  expect_identical(transcripts$score, rep(NA_real_, 2L))
  expect_identical(transcripts$exon_frames, rep(NA_character_, 2L))
  expect_identical(exons$transcript_id, rep(c("plus_tx", "minus_tx"), c(2L, 3L)))
  expect_identical(exons$gene_id, exons$transcript_id)
  expect_identical(exons$exon_number, c(1L, 2L, 3L, 2L, 1L))
  expect_identical(exons$exon_start, c(1L, 31L, 51L, 71L, 87L))
  expect_identical(exons$exon_end, c(1L, 40L, 59L, 78L, 92L))
  expect_identical(exons$exon_frame, rep(NA_integer_, 5L))
  expect_equal(nrow(annotation$validation$invalid_records), 0L)
  expect_error(read_genepred(path, format = "genePredExt", verbose = FALSE), "at least 15 columns")
})

test_that("GenePred readers preserve two-value exon lists without trailing commas", {
  path <- tempfile(fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  standard <- "two_exons\tchr1\t+\t20\t40\t20\t38\t2\t20,30\t25,40"
  rows <- c(genePred = standard,
    genePredExt = paste0(standard, "\t2.5\tgene_two\tcmpl\tcmpl\t0,2"))
  frames <- list(genePred = rep(NA_integer_, 2L), genePredExt = c(0L, 2L))
  scores <- c(genePred = NA_real_, genePredExt = 2.5)

  for (format in names(rows)) {
    writeLines(rows[[format]], path)
    annotation <- read_genepred(path, verbose = FALSE, progress = FALSE)
    exons <- as_exon_table(annotation)

    expect_identical(annotation$meta$format, format)
    expect_identical(annotation$transcripts$score, scores[[format]])
    expect_identical(exons$transcript_id, rep("two_exons", 2L))
    expect_identical(exons$exon_start, c(21L, 31L))
    expect_identical(exons$exon_end, c(25L, 40L))
    expect_identical(exons$exon_frame, frames[[format]])
    expect_equal(nrow(annotation$validation$invalid_records), 0L)
  }
})

test_that("GenePredExt preserves closed intervals and frames for valid neighboring rows", {
  path <- tempfile(fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(c(
    "iso_1\tchr1\t+\t1\t20\t4\t17\t2\t1,11,\t5,20,\t23\tgroup_a\tcmpl\tincmpl\t0,2,",
    "iso_2\tchr1\t+\t2\t26\t5\t18\t2\t2,12,\t6,26,\t17\tgroup_a\tincmpl\tcmpl\t1,",
    "ungrouped_tx\tchr2\t-\t41\t60\t41\t58\t2\t41,51,\t47,60,\t8\t\tcmpl\tcmpl\t2,1,"
  ), path)

  annotation <- read_genepred(path, coordinate = "granges", verbose = FALSE, progress = FALSE)
  transcripts <- as_transcript_table(annotation)
  exons <- as_exon_table(annotation)
  genes <- as_gene_table(annotation)

  expect_identical(annotation$meta$format, "genePredExt")
  expect_identical(annotation$meta$coordinate_input, "granges")
  expect_identical(annotation$meta$gene_col, "name2")
  expect_identical(transcripts$transcript_id, c("iso_1", "iso_2", "ungrouped_tx"))
  expect_identical(transcripts$gene_id, c("group_a", "group_a", "ungrouped_tx"))
  expect_identical(transcripts$tx_start, c(1L, 2L, 41L))
  expect_identical(transcripts$tx_end, c(20L, 26L, 60L))
  expect_identical(transcripts$cds_start, c(4L, 5L, 41L))
  expect_identical(transcripts$cds_end, c(17L, 18L, 58L))
  expect_identical(transcripts$score, c(23, 17, 8))
  expect_identical(transcripts$cds_start_stat, c("cmpl", "incmpl", "cmpl"))
  expect_identical(transcripts$cds_end_stat, c("incmpl", "cmpl", "cmpl"))
  expect_identical(exons$transcript_id, rep(transcripts$transcript_id, each = 2L))
  expect_identical(exons$exon_number, c(1L, 2L, 1L, 2L, 2L, 1L))
  expect_identical(exons$exon_start, c(1L, 11L, 2L, 12L, 41L, 51L))
  expect_identical(exons$exon_end, c(5L, 20L, 6L, 26L, 47L, 60L))
  expect_identical(exons$exon_frame, c(0L, 2L, NA_integer_, NA_integer_, 2L, 1L))
  expect_identical(genes$gene_id, c("group_a", "ungrouped_tx"))
  expect_identical(genes$n_transcripts, c(2L, 1L))
  expect_identical(genes$gene_start, c(1L, 41L))
  expect_identical(genes$gene_end, c(26L, 60L))
})

test_that("mismatched GenePred exon lists report the bad row and preserve valid rows", {
  clean_path <- tempfile(fileext = ".genePredExt")
  mixed_path <- tempfile(fileext = ".genePredExt")
  on.exit(unlink(c(clean_path, mixed_path)), add = TRUE)
  rows <- c(
    "plus_tx\tchr1\t+\t0\t40\t3\t38\t2\t0,30,\t10,40,\t20\tgene_a\tcmpl\tcmpl\t0,2,",
    "bad_tx\tchr1\t+\t50\t90\t51\t88\t1\t50,80,\t90,\t15\tgene_bad\tcmpl\tcmpl\t0,1,",
    "minus_tx\tchr2\t-\t100\t142\t102\t141\t2\t100,130,\t109,142,\t19\tgene_b\tcmpl\tcmpl\t2,1,"
  )
  writeLines(rows[c(1L, 3L)], clean_path)
  writeLines(rows, mixed_path)

  clean <- read_genepred(clean_path, verbose = FALSE, progress = FALSE)
  filtered <- read_genepred(mixed_path, verbose = FALSE, progress = FALSE)
  retained <- read_genepred(mixed_path, remove_invalid = FALSE, verbose = FALSE, progress = FALSE)
  unreported <- read_genepred(mixed_path, report_invalid = FALSE, verbose = FALSE, progress = FALSE)

  expect_identical(filtered$transcripts$transcript_id, c("plus_tx", "minus_tx"))
  for (table in c("transcripts", "exons", "genes", "data")) {
    expect_identical(as.data.frame(filtered[[table]]), as.data.frame(clean[[table]]))
    expect_identical(as.data.frame(unreported[[table]]), as.data.frame(filtered[[table]]))
  }
  expect_identical(filtered$validation$invalid_records$row_id, c(2L, 2L))
  expect_setequal(filtered$validation$invalid_records$reason,
    c("exonCount does not match exon records", "invalid exon coordinates"))
  expect_identical(retained$transcripts$transcript_id, c("plus_tx", "bad_tx", "minus_tx"))
  expect_identical(as.data.frame(retained$validation$invalid_records),
    as.data.frame(filtered$validation$invalid_records))
  expect_equal(nrow(unreported$validation$invalid_records), 0L)
})

test_that("GenePred validates mismatched exon lists with bare integer endpoints", {
  path <- tempfile(fileext = ".genePred")
  on.exit(unlink(path), add = TRUE)
  writeLines(c(
    "good_tx\tchr1\t+\t0\t10\t0\t9\t1\t0\t10",
    "bad_tx\tchr1\t+\t20\t40\t20\t38\t2\t20,30,\t40"
  ), path)

  annotation <- read_genepred(path, verbose = FALSE, progress = FALSE)
  exons <- as_exon_table(annotation)
  expect_identical(as_transcript_table(annotation)$transcript_id, "good_tx")
  expect_identical(exons$transcript_id, "good_tx")
  expect_identical(exons$exon_start, 1L)
  expect_identical(exons$exon_end, 10L)
  expect_identical(annotation$validation$invalid_records$row_id, 2L)
  expect_identical(annotation$validation$invalid_records$reason, "invalid exon coordinates")
})

test_that("GenePredExt accepts identifier overrides and rejects missing columns", {
  path <- tempfile(fileext = ".genePredExt")
  on.exit(unlink(path), add = TRUE)
  writeLines("tx_one\tchr1\t-\t0\t10\t0\t9\t1\t0\t10\t0\talternate_one\tcmpl\tcmpl\t0", path)

  annotation <- read_genepred(path, format = "genePredExt", gene_col = "name",
    transcript_col = "name2", verbose = FALSE, progress = FALSE)
  expect_identical(annotation$meta$transcript_col, "name2")
  expect_identical(as_transcript_table(annotation)$transcript_id, "alternate_one")
  expect_identical(as_transcript_table(annotation)$gene_id, "tx_one")
  expect_identical(as_exon_table(annotation)$exon_frame, 0L)
  expect_error(read_genepred(path, transcript_col = "unknown", verbose = FALSE), "transcript_col")

  writeLines("tx_one\tchr1\t-\t0\t10\t0\t9\t1\t0", path)
  expect_error(read_genepred(path, verbose = FALSE), "at least 10 columns")
})
