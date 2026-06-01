test_that("example annotation files are readable and internally consistent", {
  gp <- read_genepred(gtr_extdata("example.genePredExt"), format = "genePredExt", verbose = FALSE)
  gff <- read_gff(gtr_extdata("example_annotation.gff3"), verbose = FALSE)
  gtf <- read_gtf(gtr_extdata("example_annotation.gtf"), verbose = FALSE)
  bed <- read_bed(gtr_extdata("example_features.bed"), verbose = FALSE)

  expect_s3_class(gp, "GenePred")
  expect_s3_class(gff, "Feature")
  expect_s3_class(gtf, "Feature")
  expect_s3_class(bed, "Feature")

  expect_gte(nrow(gp$genes), 90L)
  expect_gte(nrow(gp$transcripts), 90L)
  expect_gte(nrow(gp$exons), 90L)
  expect_true("GeneA" %in% gp$genes$gene_id)
  expect_true("TxA1" %in% gp$transcripts$transcript_id)

  expect_gte(nrow(gff$genes), 90L)
  expect_gte(nrow(gtf$genes), 90L)
  expect_gte(nrow(bed$data), 90L)
})

test_that("retrieve_feature applies chromosome and range filters strictly", {
  gp <- read_genepred(gtr_extdata("example.genePredExt"), format = "genePredExt", verbose = FALSE)

  sub <- retrieve_feature(gp, chrom = "chr1", start = 1, end = 1200, as = "Feature")
  expect_s3_class(sub, "GenePred")
  expect_true(nrow(sub$genes) > 0L)
  expect_true(all(sub$genes$chrom == "chr1"))
  expect_true(all(sub$transcripts$chrom == "chr1"))
  expect_true(all(sub$exons$chrom == "chr1"))
  expect_true("GeneA" %in% sub$genes$gene_id)
})

test_that("GFF3 parser handles multi-transcript target-gene annotations", {
  gff_lines <- c(
    "##gff-version 3",
    "11\tmaker\tgene\t11046498\t11048885\t.\t+\t.\tID=GmJD_11G0141500;Name=GmJD_11G0141500",
    "11\tmaker\tmRNA\t11046498\t11048885\t.\t+\t.\tID=GmJD_11T0141500.1;Parent=GmJD_11G0141500;Name=GmJD_11T0141500.1",
    "11\tmaker\tfive_prime_UTR\t11046498\t11046725\t.\t+\t.\tID=GmJD_11T0141500.1.five_prime_UTR1;Parent=GmJD_11T0141500.1",
    "11\tmaker\texon\t11046498\t11048885\t.\t+\t.\tID=GmJD_11T0141500.1.exon1;Parent=GmJD_11T0141500.1",
    "11\tmaker\tmRNA\t11046578\t11048671\t.\t+\t.\tID=GmJD_11T0141500.2;Parent=GmJD_11G0141500;Name=GmJD_11T0141500.2",
    "11\tmaker\texon\t11046578\t11048671\t.\t+\t.\tID=GmJD_11T0141500.2.exon1;Parent=GmJD_11T0141500.2",
    "11\tmaker\tCDS\t11046579\t11048555\t.\t+\t0\tID=GmJD_11T0141500.2.CDS1;Parent=GmJD_11T0141500.2",
    "11\tmaker\tmRNA\t11047699\t11048770\t.\t+\t.\tID=GmJD_11T0141500.5;Parent=GmJD_11G0141500;Name=GmJD_11T0141500.5",
    "11\tmaker\texon\t11047699\t11048770\t.\t+\t.\tID=GmJD_11T0141500.5.exon1;Parent=GmJD_11T0141500.5",
    "11\tmaker\tthree_prime_UTR\t11048556\t11048770\t.\t+\t.\tID=GmJD_11T0141500.5.three_prime_UTR1;Parent=GmJD_11T0141500.5"
  )
  f <- tempfile(fileext = ".gff3")
  writeLines(gff_lines, f)

  x <- read_gff(f, verbose = FALSE)
  expect_s3_class(x, "Feature")
  expect_equal(nrow(x$genes), 1L)
  expect_gte(nrow(x$transcripts), 3L)
  expect_true(all(x$transcripts$gene_id == "GmJD_11G0141500"))
  expect_true(all(x$exons$gene_id == "GmJD_11G0141500"))
})

test_that("GTF parser can construct exon-like ranges from CDS-only compact annotations", {
  gtf_lines <- c(
    'chr1\tsrc\tgene\t100\t500\t.\t+\t.\tgene_id "GeneX"; gene_name "GeneX";',
    'chr1\tsrc\ttranscript\t100\t500\t.\t+\t.\tgene_id "GeneX"; transcript_id "TxX1"; gene_name "GeneX";',
    'chr1\tsrc\tCDS\t120\t220\t.\t+\t0\tgene_id "GeneX"; transcript_id "TxX1"; exon_number "1";',
    'chr1\tsrc\tCDS\t300\t450\t.\t+\t2\tgene_id "GeneX"; transcript_id "TxX1"; exon_number "2";'
  )
  f <- tempfile(fileext = ".gtf")
  writeLines(gtf_lines, f)

  x <- read_gtf(f, verbose = FALSE)
  expect_s3_class(x, "Feature")
  expect_equal(nrow(x$genes), 1L)
  expect_equal(nrow(x$transcripts), 1L)
  expect_gte(nrow(x$exons), 2L)
})

test_that("file size formatter avoids 0.00 MB for small files", {
  expect_equal(GeneTrackR:::format_file_size(523), "523 B")
  expect_match(GeneTrackR:::format_file_size(12 * 1024), "KB")
  expect_match(GeneTrackR:::format_file_size(2 * 1024^2), "MB")
})
