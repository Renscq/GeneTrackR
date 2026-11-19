test_that("feature summaries use inclusive lengths and explicit grouping", {
  object <- Feature(data.frame(
    chrom = c("chr2", "chr1", "chr1"), start = c(10L, 1L, 20L),
    end = c(14L, 10L, 39L), type = c("peak", "peak", "enhancer")
  ))
  before <- serialize(object, NULL)
  observed <- summary_feature(object)
  expect_identical(observed$chrom, c("chr1", "chr1", "chr2"))
  expect_identical(observed$type, c("enhancer", "peak", "peak"))
  expect_identical(observed$n_features, c(1L, 1L, 1L))
  expect_identical(observed$median_feature_length, c(20, 10, 5))
  fallback <- summary_feature(object, by = "not_a_column")
  expect_identical(fallback$n_features, c(2L, 1L))
  expect_identical(fallback$median_feature_length, c(15, 5))
  expect_equal(nrow(summary_feature(object, chrom = "chr9")), 0L)
  expect_identical(serialize(object, NULL), before)
  expect_error(summary_feature(list()), "Feature-compatible")
  expect_error(summary_feature(object, level = "not_a_level"), "arg")
})

test_that("gene model summaries preserve level-specific counts and empty schemas", {
  object <- coverage_gene_model()
  genes <- summary_feature(object, level = "gene")
  expect_identical(genes$strand, c("+", "-"))
  expect_identical(genes$n_genes, c(2L, 1L))
  expect_identical(genes$n_coding_genes, c(1L, 1L))
  expect_identical(genes$n_noncoding_genes, c(1L, 0L))
  expect_identical(genes$median_gene_length, c(8, 12))
  transcripts <- summary_feature(object, level = "transcript")
  expect_identical(transcripts$n_transcripts, c(2L, 1L))
  expect_identical(transcripts$median_transcript_length, c(8, 12))
  exons <- summary_feature(object, level = "exon")
  expect_identical(exons$n_exons, c(3L, 2L))
  expect_identical(exons$median_exon_length, c(4, 4.5))
  for (level in c("gene", "transcript", "exon")) {
    empty <- summary_feature(object, level = level, chrom = "chr9")
    expect_equal(nrow(empty), 0L)
    expect_true(all(c("chrom", "strand") %in% names(empty)))
  }
})

test_that("variant validation identifies independent record defects without mutation", {
  object <- VariantTrack(data.frame(
    chrom = "chr1", pos = c(10L, 20L, 30L),
    variant_id = c("v1", "v2", "v3"), ref = "A", alt = "G"
  ))
  expect_equal(nrow(validate_vcf(object)$invalid_records), 0L)
  object$data$chrom[1L] <- ""
  object$data$pos[2L] <- 0L
  object$data$ref[2L] <- NA_character_
  object$data$start[2L] <- 2L
  object$data$end[2L] <- 1L
  object$data$alt[3L] <- ""
  object$data$variant_id[3L] <- "v1"
  object$data$variant_type[3L] <- ""
  before <- serialize(object, NULL)
  observed <- validate_vcf(object)
  expect_setequal(observed$invalid_records$reason, c(
    "empty chromosome", "invalid variant position", "empty REF allele",
    "invalid start/end coordinates", "pos does not match start/end",
    "empty ALT allele", "duplicated variant_id", "empty variant_type"
  ))
  expect_identical(observed$invalid_records$row_id[
    observed$invalid_records$reason == "duplicated variant_id"
  ], 3L)
  expect_identical(observed$warnings, character())
  expect_identical(serialize(object, NULL), before)
  empty <- VariantTrack(data.frame(chrom = character(), pos = integer()))
  expect_equal(nrow(validate_vcf(empty)$invalid_records), 0L)
  expect_error(validate_vcf(list()), "VariantTrack")
})

test_that("feature validation separates invalid records from missing-type warnings", {
  object <- Feature(data.frame(chrom = "chr1", start = c(1L, 10L),
                               end = c(5L, 15L), feature_id = c("a", "b")))
  object$data$chrom[1L] <- ""
  object$data$start[1L] <- 6L
  object$data$strand[2L] <- "invalid"
  object$data$feature_id[2L] <- "a"
  object$data$type[2L] <- ""
  before <- serialize(object, NULL)
  observed <- validate_feature(object, check_gene_model = FALSE)
  expect_setequal(observed$invalid_records$reason,
                   c("empty chromosome", "invalid coordinates", "invalid strand", "duplicated feature_id"))
  expect_identical(observed$warnings, "Some feature records have missing feature type.")
  expect_identical(serialize(object, NULL), before)
  expect_error(validate_feature(list()), "Feature-compatible")
})

test_that("gene-model validation reports geometry and count defects without mutation", {
  object <- coverage_gene_model()
  object$transcripts$tx_start[1L] <- 15L
  object$transcripts$cds_start[2L] <- 20L
  object$transcripts$exon_count[3L] <- 2L
  object$transcripts$strand[3L] <- "invalid"
  object$exons$exon_start[2L] <- 4L
  object$exons$exon_end[4L] <- 34L
  object$exons$exon_start[5L] <- 45L
  before <- serialize(object, NULL)
  observed <- validate_feature(object)
  expect_true(all(c(
    "invalid transcript coordinates", "invalid CDS coordinates", "invalid strand",
    "exonCount does not match exon records", "invalid exon coordinates",
    "exon outside transcript region", "overlapping exons"
  ) %in% observed$invalid_records$reason))
  expect_identical(serialize(object, NULL), before)
})

test_that("gene-model validation warns about cross-chromosome and strand gene assignments", {
  object <- coverage_gene_model()
  object$transcripts$gene_id[seq_len(2L)] <- "shared"
  object$transcripts$chrom[2L] <- "chr2"
  observed <- validate_feature(object)
  expect_setequal(observed$warnings, c(
    "Some gene IDs are assigned to multiple chromosomes.",
    "Some gene IDs are assigned to multiple strands."
  ))
})

test_that("internal variant merging honors conflict policy and source provenance", {
  first <- VariantTrack(data.frame(chrom = c("chr2", "chr1"), pos = c(20L, 5L),
                                   variant_id = c("dup", "unique"), ref = "A", alt = "G"))
  second <- VariantTrack(data.frame(chrom = "chr1", pos = 8L, variant_id = "dup",
                                    ref = "C", alt = "T"))
  before <- serialize(list(first, second), NULL)
  merged <- GeneTrackR:::merge_vcf(first, second, source_names = c("first", "second"))
  expect_identical(merged$data$pos, c(5L, 8L, 20L))
  expect_identical(merged$data$track_source, c("first", "second", "first"))
  kept <- GeneTrackR:::merge_vcf(first, second, conflict = "keep_first")
  expect_identical(kept$data$variant_id, c("unique", "dup"))
  expect_identical(kept$data$pos, c(5L, 20L))
  renamed <- GeneTrackR:::merge_vcf(first, second, conflict = "rename",
                                   source_names = c("first", "second"))
  expect_identical(renamed$data$variant_id, c("first_unique", "second_dup", "first_dup"))
  unsorted <- GeneTrackR:::merge_vcf(first, second, sort = FALSE)
  # VariantTrack construction still applies its public canonical ordering.
  expect_identical(unsorted$data$pos, c(5L, 8L, 20L))
  expect_error(GeneTrackR:::merge_vcf(first, second, conflict = "error"), "Duplicated")
  expect_error(GeneTrackR:::merge_vcf(), "At least one")
  expect_error(GeneTrackR:::merge_vcf(first, list()), "VariantTrack")
  expect_error(GeneTrackR:::merge_vcf(first, second, source_names = "one"), "source_names")
  expect_error(GeneTrackR:::merge_vcf(first, second, rename_prefix = "one"), "rename_prefix")
  expect_identical(serialize(list(first, second), NULL), before)
})

test_that("phenotype result printing returns invisibly and renders a real figure", {
  no_figure <- structure(list(figure = NULL), class = "GeneTrackRPhenoPlot")
  empty_result <- withVisible(print(no_figure))
  expect_identical(empty_result$value, no_figure)
  expect_false(empty_result$visible)
  figure <- ggplot2::ggplot(data.frame(x = c(1, 2), y = c(3, 4)),
                            ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
  object <- structure(list(figure = figure), class = "GeneTrackRPhenoPlot")
  output <- tempfile(fileext = ".png")
  grDevices::png(output, width = 400, height = 300)
  device <- grDevices::dev.cur()
  on.exit({
    if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
    unlink(output)
  }, add = TRUE)
  observed <- withVisible(print(object))
  grDevices::dev.off(device)
  expect_identical(observed$value, object)
  expect_false(observed$visible)
  expect_identical(readBin(output, "raw", n = 8L),
                   as.raw(c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
})
