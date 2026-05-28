test_that("read_genepred reads example genePredExt", {
  file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
  gp <- read_genepred(file, format = "genePredExt", verbose = FALSE)

  expect_s3_class(gp, "GenePred")
  expect_equal(nrow(gp$genes), 4L)
  expect_equal(nrow(gp$transcripts), 5L)
  expect_equal(nrow(gp$exons), 10L)
  expect_true(all(c("GeneA", "GeneB", "GeneC", "GeneD") %in% gp$genes$gene_id))
})

test_that("slice_genepred supports region slicing", {
  file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
  gp <- read_genepred(file, format = "genePredExt", verbose = FALSE)
  sub <- slice_genepred(gp, chrom = "chr1", start = 1, end = 1000, mode = "overlap")

  expect_s3_class(sub, "GenePred")
  expect_equal(unique(sub$genes$gene_id), "GeneA")
})

test_that("summary_genepred returns stable summary tables", {
  file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
  gp <- read_genepred(file, format = "genePredExt", verbose = FALSE)

  gene_sum <- summary_genepred(gp, level = "gene")
  tx_sum <- summary_genepred(gp, level = "transcript")
  exon_sum <- summary_genepred(gp, level = "exon")

  expect_true(nrow(gene_sum) > 0L)
  expect_true(nrow(tx_sum) > 0L)
  expect_true(nrow(exon_sum) > 0L)
})


test_that("search_gene and search_transcript find annotation records", {
  file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
  gp <- read_genepred(file, format = "genePredExt", verbose = FALSE)

  genes <- search_gene(gp, "GeneA")
  tx <- search_transcript(gp, "TxA")

  expect_equal(genes$gene_id, "GeneA")
  expect_true(all(c("TxA1", "TxA2") %in% tx$transcript_id))
})
