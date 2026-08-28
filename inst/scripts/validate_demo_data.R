# Author: Rensc
# Date: 2026-08-23
# Version: dev009
# Function: Validate structural invariants and designed truth of the GeneTrackR demo dataset
# Input: inst/scripts/demo_model/*.tsv and inst/extdata/gtr_demo_* files
# Output: Validation messages; stops on any failed invariant

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

read_tsv <- function(file) {
  utils::read.delim(file, sep = "	", quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
}

expect_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

parse_int_list <- function(x) {
  as.integer(strsplit(as.character(x), ",", fixed = TRUE)[[1L]])
}

get_cds_positions <- function(tx_row) {
  if (is.na(tx_row$cds_start) || is.na(tx_row$cds_end)) {
    return(integer())
  }
  starts <- parse_int_list(tx_row$exon_starts)
  ends <- parse_int_list(tx_row$exon_ends)
  cds_start <- as.integer(tx_row$cds_start)
  cds_end <- as.integer(tx_row$cds_end)
  segments <- lapply(seq_along(starts), function(i) {
    start <- max(starts[i], cds_start)
    end <- min(ends[i], cds_end)
    if (start > end) return(NULL)
    c(start, end)
  })
  segments <- Filter(Negate(is.null), segments)
  if (length(segments) == 0L) return(integer())
  if (tx_row$strand == "+") {
    segment_starts <- vapply(segments, function(x) x[1L], integer(1L))
    segments <- segments[order(segment_starts)]
    return(unlist(lapply(segments, function(x) seq.int(x[1L], x[2L])), use.names = FALSE))
  }
  segment_starts <- vapply(segments, function(x) x[1L], integer(1L))
  segments <- segments[order(segment_starts, decreasing = TRUE)]
  unlist(lapply(segments, function(x) seq.int(x[2L], x[1L])), use.names = FALSE)
}

script_dir <- get_script_dir()
model_dir <- file.path(script_dir, "demo_model")
ext_dir <- normalizePath(file.path(script_dir, "..", "extdata"), mustWork = FALSE)

required_files <- c(
  "gtr_demo.genePredExt", "gtr_demo.gff3", "gtr_demo.gtf", "gtr_demo_features.bed",
  "gtr_demo_variants.vcf", "gtr_demo_pheno.tsv",
  "gtr_demo_rnaseq_plus.bedgraph", "gtr_demo_rnaseq_minus.bedgraph",
  "gtr_demo_riboseq_plus.bedgraph", "gtr_demo_riboseq_minus.bedgraph",
  "gtr_demo.chrom.sizes"
)
expect_true(all(file.exists(file.path(ext_dir, required_files))), "One or more generated demo files are missing.")
legacy_example_files <- c(
  "example.genePredExt",
  "example_annotation.gff3",
  "example_annotation.gtf",
  "example_features.bed",
  "example_haplotype.vcf",
  "example_pheno.tsv",
  "example_signal_A.bedgraph",
  "example_signal_B.bedgraph",
  "example_variants.vcf",
  "example_variants_NC12.vcf"
)
expect_true(!any(file.exists(file.path(ext_dir, legacy_example_files))), "Legacy example files should not remain in inst/extdata.")

chromosomes <- read_tsv(file.path(model_dir, "chromosomes.tsv"))
genes <- read_tsv(file.path(model_dir, "genes.tsv"))
transcripts <- read_tsv(file.path(model_dir, "transcripts.tsv"))
samples <- read_tsv(file.path(model_dir, "samples.tsv"))
variants <- read_tsv(file.path(model_dir, "variants.tsv"))
signal_design <- read_tsv(file.path(model_dir, "signal_design.tsv"))
pheno <- read_tsv(file.path(ext_dir, "gtr_demo_pheno.tsv"))
rnaseq_plus <- utils::read.delim(
  file.path(ext_dir, "gtr_demo_rnaseq_plus.bedgraph"),
  header = FALSE, sep = "	", stringsAsFactors = FALSE,
  col.names = c("chrom", "start", "end", "value")
)
rnaseq_minus <- utils::read.delim(
  file.path(ext_dir, "gtr_demo_rnaseq_minus.bedgraph"),
  header = FALSE, sep = "	", stringsAsFactors = FALSE,
  col.names = c("chrom", "start", "end", "value")
)
riboseq_plus <- utils::read.delim(
  file.path(ext_dir, "gtr_demo_riboseq_plus.bedgraph"),
  header = FALSE, sep = "	", stringsAsFactors = FALSE,
  col.names = c("chrom", "start", "end", "value")
)
riboseq_minus <- utils::read.delim(
  file.path(ext_dir, "gtr_demo_riboseq_minus.bedgraph"),
  header = FALSE, sep = "	", stringsAsFactors = FALSE,
  col.names = c("chrom", "start", "end", "value")
)

rnaseq <- rbind(rnaseq_plus, rnaseq_minus)
riboseq <- rbind(riboseq_plus, riboseq_minus)

expect_true(nrow(chromosomes) == 2L, "Expected exactly 2 chromosomes.")
expect_true(nrow(genes) == 20L, "Expected exactly 20 genes.")
expect_true(nrow(transcripts) == 24L, "Expected exactly 24 transcripts.")
expect_true(nrow(samples) == 36L, "Expected exactly 36 samples.")
expect_true(nrow(variants) == 56L, "Expected exactly 56 variants.")
expect_true(nrow(signal_design) == 24L, "Expected one signal-design row per transcript.")
expect_true(nrow(pheno) == 36L, "Expected exactly 36 phenotype rows.")
expect_true(setequal(signal_design$transcript_id, transcripts$transcript_id), "Signal-design transcript IDs do not match the transcript model.")
signal_tx <- merge(
  signal_design,
  transcripts[, c("transcript_id", "gene_id", "gene_type")],
  by = "transcript_id",
  suffixes = c("_signal", "_model"),
  sort = FALSE
)
expect_true(all(signal_tx$gene_id_signal == signal_tx$gene_id_model), "Signal-design gene IDs do not match transcript gene IDs.")
expect_true(all(signal_design$riboseq_weight[!as.logical(signal_design$is_primary)] == 0), "Alternative demo isoforms must not contribute Ribo-seq density.")
expect_true(all(signal_tx$riboseq_weight[signal_tx$gene_type == "lncRNA"] == 0), "lncRNA transcripts must have zero Ribo-seq weight.")

coding_tx <- transcripts[transcripts$gene_type == "protein_coding", , drop = FALSE]
coding_cds_lengths <- vapply(seq_len(nrow(coding_tx)), function(i) {
  length(get_cds_positions(coding_tx[i, , drop = FALSE]))
}, integer(1L))
expect_true(
  all(coding_cds_lengths > 0L & coding_cds_lengths %% 3L == 0L),
  "Every protein-coding demo transcript must have a non-zero CDS length divisible by 3."
)
expect_true(setequal(pheno$sample_id, samples$sample_id), "Phenotype and VCF-design sample IDs do not match.")
expect_true(!identical(pheno$sample_id, samples$sample_id), "Phenotype rows should intentionally differ from VCF sample order.")
expect_true(!anyDuplicated(samples$sample_id), "Duplicated sample IDs were found in sample design.")
expect_true(!anyDuplicated(variants$variant_id), "Duplicated variant IDs were found in variant design.")

hap_counts <- table(samples$core_haplotype)
expect_true(length(hap_counts) == 4L && all(hap_counts == 9L), "Expected four balanced core haplotypes with 9 samples each.")

ld_ids <- paste0("varLD", sprintf("%02d", 1:12))
ld <- variants[variants$variant_id %in% ld_ids, , drop = FALSE]
expect_true(nrow(ld) == 12L, "Expected exactly 12 variants in the primary LD-gradient example.")
expect_true(all(ld$pattern[1:6] == "p13"), "varLD01-varLD06 must retain the perfect p13 LD core.")
expect_true(identical(ld$pattern[7:12], c("ldgrad02", "ldgrad04", "ldgrad06", "ldgrad09", "ldgrad12", "ldgrad18")), "The downstream LD-gradient patterns are not in the designed order.")
expect_true(all(ld$pos[7:12] > 12352000L & ld$pos[7:12] < 12356001L), "The added LD-gradient variants must remain outside the GeneA and GeneB gene bodies.")

pair <- variants[variants$chrom == "chr2" & variants$pos >= 16995001L & variants$pos <= 17006000L, , drop = FALSE]
expect_true(identical(pair$variant_id, c("varPair01", "varPair02")), "GeneT two-variant LD region must contain exactly varPair01 and varPair02.")

expect_true(all(c("seed_weight", "protein_content", "plant_height", "flowering_time", "flower_color") %in% names(pheno)), "Phenotype table is missing designed traits.")


varA03_design <- variants[variants$variant_id == "varA03", , drop = FALSE]
expect_true(nrow(varA03_design) == 1L && identical(varA03_design$pattern, "p13"), "varA03 must follow the p13 genotype-cluster split.")

pheno_design <- merge(samples, pheno, by = "sample_id", sort = FALSE)
seed_means <- tapply(pheno_design$seed_weight, pheno_design$hap_group, mean)
protein_means <- tapply(pheno_design$protein_content, pheno_design$hap_group, mean)
height_means <- tapply(pheno_design$plant_height, pheno_design$hap_group, mean)
expect_true(isTRUE(all.equal(as.numeric(seed_means), c(20, 22, 30, 28), tolerance = 1e-8)), "seed_weight group means do not match the hierarchical haplotype design.")
expect_true(isTRUE(all.equal(as.numeric(protein_means), c(38, 38, 44, 44), tolerance = 1e-8)), "protein_content group means do not match the genotype-nearest refinement design.")
expect_true(isTRUE(all.equal(as.numeric(height_means), c(100, 102, 108, 106), tolerance = 1e-8)), "plant_height group means do not match the hierarchical haplotype design.")


# The negative-control trait uses the same within-haplotype pattern in all four groups.
flower_patterns <- split(pheno_design$flowering_time, pheno_design$hap_group)
flower_patterns <- lapply(flower_patterns, sort)
expect_true(all(vapply(flower_patterns[-1L], identical, logical(1L), flower_patterns[[1L]])), "flowering_time negative-control distributions differ across haplotypes.")

# Signal-track invariants.
expect_true(all(rnaseq_plus$start < rnaseq_plus$end), "RNA-seq plus bedGraph contains invalid intervals.")
expect_true(all(rnaseq_minus$start < rnaseq_minus$end), "RNA-seq minus bedGraph contains invalid intervals.")
expect_true(all(rnaseq_plus$value > 0), "RNA-seq plus bedGraph must contain positive coverage values only.")
expect_true(all(rnaseq_minus$value > 0), "RNA-seq minus bedGraph must contain positive coverage values only.")
expect_true(all(riboseq_plus$end - riboseq_plus$start == 1L), "Every Ribo-seq plus bedGraph interval must represent exactly one genomic base.")
expect_true(all(riboseq_minus$end - riboseq_minus$start == 1L), "Every Ribo-seq minus bedGraph interval must represent exactly one genomic base.")
expect_true(all(riboseq_plus$value > 0), "Ribo-seq plus bedGraph must contain positive counts only.")
expect_true(all(riboseq_minus$value > 0), "Ribo-seq minus bedGraph must contain positive counts only.")
expect_true(all(riboseq_plus$value == round(riboseq_plus$value)), "Ribo-seq plus bedGraph must contain integer P-site counts.")
expect_true(all(riboseq_minus$value == round(riboseq_minus$value)), "Ribo-seq minus bedGraph must contain integer P-site counts.")
expect_true(!anyDuplicated(paste(riboseq_plus$chrom, riboseq_plus$start, sep = ":")), "Ribo-seq plus bedGraph contains duplicated genomic positions.")
expect_true(!anyDuplicated(paste(riboseq_minus$chrom, riboseq_minus$start, sep = ":")), "Ribo-seq minus bedGraph contains duplicated genomic positions.")

# Strand-specific signal checks.
expect_true(any(rnaseq_plus$chrom == "chr1" & rnaseq_plus$start < 12341000L & rnaseq_plus$end > 12340000L), "RNA-seq plus signal is missing from the plus-strand GeneA exon region.")
expect_true(!any(rnaseq_minus$chrom == "chr1" & rnaseq_minus$start < 12341000L & rnaseq_minus$end > 12340000L), "RNA-seq minus signal should be absent from the plus-strand GeneA exon region.")
expect_true(any(rnaseq_minus$chrom == "chr1" & rnaseq_minus$start < 12357500L & rnaseq_minus$end > 12356000L), "RNA-seq minus signal is missing from the minus-strand GeneB exon region.")
expect_true(!any(rnaseq_plus$chrom == "chr1" & rnaseq_plus$start < 12357500L & rnaseq_plus$end > 12356000L), "RNA-seq plus signal should be absent from the minus-strand GeneB exon region.")
expect_true(any(riboseq_plus$chrom == "chr1" & riboseq_plus$start >= 12340500L & riboseq_plus$end <= 12351500L), "Ribo-seq plus signal is missing from the plus-strand GeneA CDS region.")
expect_true(!any(riboseq_minus$chrom == "chr1" & riboseq_minus$start >= 12340500L & riboseq_minus$end <= 12351500L), "Ribo-seq minus signal should be absent from the plus-strand GeneA CDS region.")
expect_true(any(riboseq_minus$chrom == "chr1" & riboseq_minus$start >= 12356500L & riboseq_minus$end <= 12366000L), "Ribo-seq minus signal is missing from the minus-strand GeneB CDS region.")
expect_true(!any(riboseq_plus$chrom == "chr1" & riboseq_plus$start >= 12356500L & riboseq_plus$end <= 12366000L), "Ribo-seq plus signal should be absent from the minus-strand GeneB CDS region.")

# GeneA exon 1 has RNA-seq coverage, whereas the following intron is signal-free.
genea_exon_rna <- rnaseq$chrom == "chr1" & rnaseq$start < 12341000L & rnaseq$end > 12340000L
genea_intron_rna <- rnaseq$chrom == "chr1" & rnaseq$start < 12342000L & rnaseq$end > 12341000L
expect_true(any(genea_exon_rna), "RNA-seq signal is missing from the GeneA exon region.")
expect_true(!any(genea_intron_rna), "RNA-seq signal should be absent from the designed GeneA intron.")

# The non-coding GeneI is transcribed in RNA-seq but must have no Ribo-seq CDS signal.
genei_rna <- rnaseq$chrom == "chr1" & rnaseq$start < 15006000L & rnaseq$end > 15000000L
genei_ribo <- riboseq$chrom == "chr1" & riboseq$start < 15006000L & riboseq$end > 15000000L
expect_true(any(genei_rna), "RNA-seq signal is missing from the designed lncRNA GeneI.")
expect_true(!any(genei_ribo), "Ribo-seq signal must not be generated for the lncRNA GeneI.")

# GeneA Ribo-seq is moderately dense P-site-like signal with an 80:10:10 frame-count ratio.
genea_tx <- transcripts[transcripts$transcript_id == "TxA1", , drop = FALSE]
genea_positions <- get_cds_positions(genea_tx)
expect_true(length(genea_positions) > 0L && length(genea_positions) %% 3L == 0L, "GeneA CDS must be a non-zero multiple of 3.")

genea_ribo <- riboseq_plus[riboseq_plus$chrom == "chr1", , drop = FALSE]
genea_ribo$pos <- genea_ribo$start + 1L
genea_counts <- numeric(length(genea_positions))
match_index <- match(genea_ribo$pos, genea_positions)
keep_match <- !is.na(match_index)
genea_counts[match_index[keep_match]] <- genea_ribo$value[keep_match]

phase <- (seq_along(genea_positions) - 1L) %% 3L
phase_totals <- vapply(0:2, function(x) sum(genea_counts[phase == x]), numeric(1L))
phase0_fraction <- phase_totals[1L] / sum(phase_totals)
occupancy <- mean(genea_counts > 0)
frame0_indices <- which(phase == 0L)
start_peak <- genea_counts[frame0_indices[1L]]
stop_peak <- genea_counts[frame0_indices[length(frame0_indices)]]
internal_indices <- frame0_indices[frame0_indices > 9L & frame0_indices < length(genea_counts) - 8L]
internal_mean <- mean(genea_counts[internal_indices][genea_counts[internal_indices] > 0])
phase_fractions <- phase_totals / sum(phase_totals)

frame_values <- lapply(0:2, function(x) genea_counts[phase == x & genea_counts > 0])
frame_unique_n <- vapply(frame_values, function(x) length(unique(x)), integer(1L))
frame_cv <- vapply(frame_values, function(x) stats::sd(x) / mean(x), numeric(1L))

expect_true(occupancy > 0.45 && occupancy < 0.65, "GeneA Ribo-seq signal should have moderate base-level occupancy.")
expect_true(all(frame_unique_n >= c(8L, 4L, 4L)), "GeneA Ribo-seq frames should contain heterogeneous non-zero count heights.")
expect_true(all(frame_cv > 0.20), "GeneA Ribo-seq frame heights should show visible within-frame variation.")
expect_true(abs(phase_fractions[1L] - 0.80) < 0.02, "GeneA Ribo-seq frame-0 count fraction should be approximately 80%.")
expect_true(abs(phase_fractions[2L] - 0.10) < 0.02, "GeneA Ribo-seq frame-1 count fraction should be approximately 10%.")
expect_true(abs(phase_fractions[3L] - 0.10) < 0.02, "GeneA Ribo-seq frame-2 count fraction should be approximately 10%.")
expect_true(abs(start_peak / internal_mean - 2) < 0.35, "GeneA Ribo-seq initiation count should be close to two times the internal frame-0 mean.")
expect_true(abs(stop_peak / internal_mean - 2) < 0.35, "GeneA Ribo-seq termination count should be close to two times the internal frame-0 mean.")

# Apply the frame-ratio and boundary-peak checks to every translated primary transcript.
translated_design <- signal_design[signal_design$riboseq_weight > 0, , drop = FALSE]
for (i in seq_len(nrow(translated_design))) {
  tx_id <- translated_design$transcript_id[i]
  tx <- transcripts[transcripts$transcript_id == tx_id, , drop = FALSE]
  tx_positions <- get_cds_positions(tx)
  expect_true(length(tx_positions) > 0L && length(tx_positions) %% 3L == 0L, paste0(tx_id, " CDS is not a valid 3n coding sequence."))

  source <- if (tx$strand[1L] == "+") riboseq_plus else riboseq_minus
  source <- source[source$chrom == tx$chrom[1L], , drop = FALSE]
  source$pos <- source$start + 1L
  tx_counts <- numeric(length(tx_positions))
  idx <- match(source$pos, tx_positions)
  matched <- !is.na(idx)
  tx_counts[idx[matched]] <- source$value[matched]

  tx_phase <- (seq_along(tx_positions) - 1L) %% 3L
  totals <- vapply(0:2, function(x) sum(tx_counts[tx_phase == x]), numeric(1L))
  fractions <- totals / sum(totals)
  expect_true(abs(fractions[1L] - 0.80) < 0.03, paste0(tx_id, " frame-0 fraction is not approximately 80%."))
  expect_true(abs(fractions[2L] - 0.10) < 0.03, paste0(tx_id, " frame-1 fraction is not approximately 10%."))
  expect_true(abs(fractions[3L] - 0.10) < 0.03, paste0(tx_id, " frame-2 fraction is not approximately 10%."))

  frame_values_tx <- lapply(0:2, function(x) tx_counts[tx_phase == x & tx_counts > 0])
  expect_true(length(unique(frame_values_tx[[1L]])) >= 5L, paste0(tx_id, " frame-0 counts are insufficiently heterogeneous."))
  expect_true(length(unique(frame_values_tx[[2L]])) >= 2L, paste0(tx_id, " frame-1 counts are insufficiently heterogeneous."))
  expect_true(length(unique(frame_values_tx[[3L]])) >= 2L, paste0(tx_id, " frame-2 counts are insufficiently heterogeneous."))

  frame0_idx <- which(tx_phase == 0L)
  internal <- tx_counts[frame0_idx[-c(1L, length(frame0_idx))]]
  internal <- internal[internal > 0]
  internal_mean_tx <- mean(internal)
  expect_true(abs(tx_counts[frame0_idx[1L]] / internal_mean_tx - 2) < 0.40, paste0(tx_id, " initiation count is not approximately 2x the internal frame-0 mean."))
  expect_true(abs(tx_counts[frame0_idx[length(frame0_idx)]] / internal_mean_tx - 2) < 0.40, paste0(tx_id, " termination count is not approximately 2x the internal frame-0 mean."))
}

vcf_lines <- readLines(file.path(ext_dir, "gtr_demo_variants.vcf"), warn = FALSE)
expect_true(any(vcf_lines == "##source=GeneTrackR_demo_v0.5.5"), "VCF source metadata does not match demo model v0.5.5.")
vcf_header <- grep("^#CHROM", vcf_lines, value = TRUE)
expect_true(length(vcf_header) == 1L, "VCF #CHROM header was not found.")
vcf_samples <- strsplit(vcf_header, "	", fixed = TRUE)[[1L]][10:45]
expect_true(identical(vcf_samples, samples$sample_id), "VCF sample order does not match the canonical sample design.")

message("[GeneTrackR demo] Validation passed.")
message("[GeneTrackR demo] chromosomes=2; genes=20; transcripts=24; samples=36; variants=50; core_haplotypes=4; signal_tracks=4.")
