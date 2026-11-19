coverage_gene_model <- function() {
  transcripts <- data.frame(
    transcript_id = c("plus_tx", "minus_tx", "noncoding_tx"),
    gene_id = c("plus_gene", "minus_gene", "noncoding_gene"),
    chrom = "chr1", strand = c("+", "-", "+"),
    tx_start = c(1L, 21L, 41L), tx_end = c(12L, 32L, 44L),
    cds_start = c(3L, 23L, 41L), cds_end = c(10L, 30L, 40L),
    gene_type = c("coding", "coding", "non-coding"),
    exon_count = c(2L, 2L, 1L)
  )
  exons <- data.frame(
    transcript_id = c("plus_tx", "plus_tx", "minus_tx", "minus_tx", "noncoding_tx"),
    gene_id = c("plus_gene", "plus_gene", "minus_gene", "minus_gene", "noncoding_gene"),
    chrom = "chr1", strand = c("+", "+", "-", "-", "+"),
    exon_number = c(1L, 2L, 2L, 1L, 1L),
    exon_start = c(1L, 9L, 21L, 29L, 41L),
    exon_end = c(5L, 12L, 25L, 32L, 44L), exon_frame = NA_integer_
  )
  GenePred(transcripts, exons)
}
