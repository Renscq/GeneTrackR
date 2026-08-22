# Author: Rensc
# Date: 2026-08-23
# Version: dev003
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
expect_true(nrow(variants) == 50L, "Expected exactly 50 variants.")
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
expect_true(setequal(pheno$sample_id, samples$sample_id), "Phenotype and VCF-design sample IDs do not match.")
expect_true(!identical(pheno$sample_id, samples$sample_id), "Phenotype rows should intentionally differ from VCF sample order.")
expect_true(!anyDuplicated(samples$sample_id), "Duplicated sample IDs were found in sample design.")
expect_true(!anyDuplicated(variants$variant_id), "Duplicated variant IDs were found in variant design.")

hap_counts <- table(samples$core_haplotype)
expect_true(length(hap_counts) == 4L && all(hap_counts == 9L), "Expected four balanced core haplotypes with 9 samples each.")

ld_ids <- paste0("varLD", sprintf("%02d", 1:6))
ld <- variants[variants$variant_id %in% ld_ids, , drop = FALSE]
expect_true(nrow(ld) == 6L && length(unique(ld$pattern)) == 1L && unique(ld$pattern) == "p13", "High-LD variants must share the p13 genotype pattern.")

pair <- variants[variants$chrom == "chr2" & variants$pos >= 16995001L & variants$pos <= 17006000L, , drop = FALSE]
expect_true(identical(pair$variant_id, c("varPair01", "varPair02")), "GeneT two-variant LD region must contain exactly varPair01 and varPair02.")

expect_true(all(c("seed_weight", "protein_content", "plant_height", "flowering_time", "flower_color") %in% names(pheno)), "Phenotype table is missing designed traits.")

# The negative-control trait uses the same within-haplotype pattern in all four groups.
pheno_design <- merge(samples, pheno, by = "sample_id", sort = FALSE)
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
expect_true(all(riboseq_plus$value > 0), "Ribo-seq plus bedGraph must contain positive density values only.")
expect_true(all(riboseq_minus$value > 0), "Ribo-seq minus bedGraph must contain positive density values only.")
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

# GeneA Ribo-seq density is anchored to the CDS start and displays a 3-nt periodic pattern away from initiation/termination peaks.
genea_ribo <- riboseq_plus[riboseq_plus$chrom == "chr1" & riboseq_plus$start >= 12340500L & riboseq_plus$end <= 12351500L, , drop = FALSE]
expect_true(nrow(genea_ribo) > 0L, "Ribo-seq signal is missing from GeneA CDS.")
genea_ribo$pos <- genea_ribo$start + 1L
start_peak <- genea_ribo$value[genea_ribo$pos == 12340501L]
stop_peak <- max(genea_ribo$value[genea_ribo$pos >= 12351498L & genea_ribo$pos <= 12351500L])
periodic <- genea_ribo[genea_ribo$pos >= 12340531L & genea_ribo$pos <= 12340990L, , drop = FALSE]
periodic$phase <- (periodic$pos - 12340501L) %% 3L
phase_means <- tapply(periodic$value, periodic$phase, mean)
expect_true(length(start_peak) == 1L && start_peak > max(periodic$value), "GeneA Ribo-seq initiation peak is not higher than the internal CDS density.")
expect_true(is.finite(stop_peak) && stop_peak > max(periodic$value), "GeneA Ribo-seq termination peak is not higher than the internal CDS density.")
expect_true(phase_means[["0"]] > phase_means[["1"]] && phase_means[["1"]] > phase_means[["2"]], "GeneA Ribo-seq signal does not show the designed 3-nt periodicity.")

vcf_header <- grep("^#CHROM", readLines(file.path(ext_dir, "gtr_demo_variants.vcf"), warn = FALSE), value = TRUE)
expect_true(length(vcf_header) == 1L, "VCF #CHROM header was not found.")
vcf_samples <- strsplit(vcf_header, "	", fixed = TRUE)[[1L]][10:45]
expect_true(identical(vcf_samples, samples$sample_id), "VCF sample order does not match the canonical sample design.")

message("[GeneTrackR demo] Validation passed.")
message("[GeneTrackR demo] chromosomes=2; genes=20; transcripts=24; samples=36; variants=50; core_haplotypes=4; signal_tracks=4.")
