# Author: Rensc
# Date: 2026-08-19
# Version: dev001
# Function: Build haplotype tables from VCF variants
# Input: VariantTrack objects, VCF files, and genomic locators
# Output: HapVariant objects

#' Build haplotypes from variants in a gene or transcript region
#'
#' @description
#' Extracts variants for a specified gene or transcript, optionally including
#' upstream and downstream flanking regions, and converts sample genotype
#' profiles into haplotypes.
#'
#' @param vcf A VariantTrack object or VCF file path.
#' @param annotation A gene annotation object used to locate `gene_id` or `transcript_id`.
#' @param gene_id Optional gene ID. Use exactly one of `gene_id` or `transcript_id`.
#' @param transcript_id Optional transcript ID. Use exactly one of `gene_id` or `transcript_id`.
#' @param upstream Upstream flanking length in bp. When `strand_aware = TRUE`, upstream is interpreted relative to the gene/transcript strand.
#' @param downstream Downstream flanking length in bp. When `strand_aware = TRUE`, downstream is interpreted relative to the gene/transcript strand.
#' @param strand_aware Logical. Whether upstream/downstream should follow strand direction. Default TRUE.
#' @param samples Optional sample names to keep.
#' @param variant_type Optional variant types to keep.
#' @param genotype_mode Genotype representation. `code` converts genotypes to compact 0/1 states, where 0 means reference genotype and 1 means any alternate allele is present. `string` converts genotypes to a single allele label; long InDel alleles are compressed as `iN`, where `N` is the displayed REF or ALT allele length, so the haplotype cell uses the same label as the corresponding REF/ALT row.
#' @param missing_genotype Missing genotype label. Default is `NA_character_`, which is displayed as `NA` in haplotype tables.
#' @param min_variant_number Minimum number of non-missing variants required for a sample. If NULL, only samples with complete non-missing genotypes across all retained variants are kept.
#' @return A HapVariant object.
#' @examples
#' vcf_file <- system.file("extdata", "example_haplotype.vcf", package = "GeneTrackR")
#' anno_file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file)
#' anno <- read_genepred(anno_file, format = "genePredExt", verbose = FALSE)
#' hap <- hap_gene_variant(vcf, annotation = anno, gene_id = "GeneA", upstream = 500, downstream = 300)
#' hap
#' hap$haplotypes
#' hap_tx <- hap_gene_variant(vcf, annotation = anno, transcript_id = "TxA1", genotype_mode = "string")
#' hap_tx$haplotypes
#' @export
hap_gene_variant <- function(vcf,
                             annotation,
                             gene_id = NULL,
                             transcript_id = NULL,
                             upstream = 0L,
                             downstream = 0L,
                             strand_aware = TRUE,
                             samples = NULL,
                             variant_type = NULL,
                             genotype_mode = c("code", "string"),
                             missing_genotype = NA_character_,
                             min_variant_number = NULL) {
  genotype_mode <- match.arg(genotype_mode)
  region <- resolve_haplotype_gene_region(
    annotation = annotation,
    gene_id = gene_id,
    transcript_id = transcript_id,
    upstream = upstream,
    downstream = downstream,
    strand_aware = strand_aware
  )

  build_haplotype_from_region(
    vcf = vcf,
    region = region,
    samples = samples,
    variant_type = variant_type,
    genotype_mode = genotype_mode,
    missing_genotype = missing_genotype,
    min_variant_number = min_variant_number
  )
}

#' Build haplotypes from variants in a genomic region
#'
#' @description
#' Extracts variants from a user-defined genomic interval and converts sample
#' genotype profiles into haplotypes.
#'
#' @param vcf A VariantTrack object or VCF file path.
#' @param chrom Chromosome name.
#' @param start Region start in 1-based closed coordinates.
#' @param end Region end in 1-based closed coordinates.
#' @param samples Optional sample names to keep.
#' @param variant_type Optional variant types to keep.
#' @param genotype_mode Genotype representation. `code` converts genotypes to compact 0/1 states, where 0 means reference genotype and 1 means any alternate allele is present. `string` converts genotypes to a single allele label; long InDel alleles are compressed as `iN`, where `N` is the displayed REF or ALT allele length, so the haplotype cell uses the same label as the corresponding REF/ALT row.
#' @param missing_genotype Missing genotype label. Default is `NA_character_`, which is displayed as `NA` in haplotype tables.
#' @param min_variant_number Minimum number of non-missing variants required for a sample. If NULL, only samples with complete non-missing genotypes across all retained variants are kept.
#' @return A HapVariant object.
#' @examples
#' vcf_file <- system.file("extdata", "example_haplotype.vcf", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file)
#' hap <- hap_region_variant(vcf, chrom = "chr1", start = 1000, end = 12000, genotype_mode = "code")
#' hap
#' hap$haplotypes
#' @export
hap_region_variant <- function(vcf,
                               chrom,
                               start,
                               end,
                               samples = NULL,
                               variant_type = NULL,
                               genotype_mode = c("code", "string"),
                               missing_genotype = NA_character_,
                               min_variant_number = NULL) {
  genotype_mode <- match.arg(genotype_mode)
  check_region(chrom, start, end)
  region <- list(
    locator = "region",
    id = paste0(as.character(chrom)[1L], ":", as.integer(start)[1L], "-", as.integer(end)[1L]),
    chrom = as.character(chrom)[1L],
    start = as.integer(start)[1L],
    end = as.integer(end)[1L],
    core_start = as.integer(start)[1L],
    core_end = as.integer(end)[1L],
    upstream = 0L,
    downstream = 0L,
    strand = NA_character_
  )

  build_haplotype_from_region(
    vcf = vcf,
    region = region,
    samples = samples,
    variant_type = variant_type,
    genotype_mode = genotype_mode,
    missing_genotype = missing_genotype,
    min_variant_number = min_variant_number
  )
}

#' Build haplotypes from variants in a gene, transcript, or genomic region
#'
#' @description
#' Compatibility wrapper for haplotype construction. For new code, prefer
#' `hap_gene_variant()` for gene/transcript queries and `hap_region_variant()`
#' for direct genomic interval queries.
#'
#' @param vcf A VariantTrack object or VCF file path.
#' @param annotation Optional gene annotation object used for `gene_id` or `transcript_id` queries.
#' @param gene_id Optional gene ID.
#' @param transcript_id Optional transcript ID.
#' @param chrom Optional chromosome name.
#' @param start Optional region start in 1-based closed coordinates.
#' @param end Optional region end in 1-based closed coordinates.
#' @param upstream Upstream flanking length in bp for gene/transcript queries.
#' @param downstream Downstream flanking length in bp for gene/transcript queries.
#' @param strand_aware Logical. Whether upstream/downstream should follow strand direction. Default TRUE.
#' @param samples Optional sample names to keep.
#' @param variant_type Optional variant types to keep.
#' @param genotype_mode Genotype representation. `code` converts genotypes to compact 0/1 states, where 0 means reference genotype and 1 means any alternate allele is present. `string` converts genotypes to a single allele label; long InDel alleles are compressed as `iN`, where `N` is the displayed REF or ALT allele length, so the haplotype cell uses the same label as the corresponding REF/ALT row.
#' @param missing_genotype Missing genotype label. Default is `NA_character_`, which is displayed as `NA` in haplotype tables.
#' @param min_variant_number Minimum number of non-missing variants required for a sample. If NULL, only samples with complete non-missing genotypes across all retained variants are kept.
#' @return A HapVariant object.
#' @examples
#' vcf_file <- system.file("extdata", "example_haplotype.vcf", package = "GeneTrackR")
#' anno_file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file)
#' anno <- read_genepred(anno_file, format = "genePredExt", verbose = FALSE)
#' hap <- hap_variant(vcf, annotation = anno, gene_id = "GeneA", genotype_mode = "code")
#' hap
#' hap_region <- hap_variant(vcf, chrom = "chr1", start = 1000, end = 12000)
#' hap_region$haplotypes
#' @export
hap_variant <- function(vcf,
                        annotation = NULL,
                        gene_id = NULL,
                        transcript_id = NULL,
                        chrom = NULL,
                        start = NULL,
                        end = NULL,
                        upstream = 0L,
                        downstream = 0L,
                        strand_aware = TRUE,
                        samples = NULL,
                        variant_type = NULL,
                        genotype_mode = c("code", "string"),
                        missing_genotype = NA_character_,
                        min_variant_number = NULL) {
  genotype_mode <- match.arg(genotype_mode)

  if (!is.null(gene_id) || !is.null(transcript_id)) {
    return(hap_gene_variant(
      vcf = vcf,
      annotation = annotation,
      gene_id = gene_id,
      transcript_id = transcript_id,
      upstream = upstream,
      downstream = downstream,
      strand_aware = strand_aware,
      samples = samples,
      variant_type = variant_type,
      genotype_mode = genotype_mode,
      missing_genotype = missing_genotype,
      min_variant_number = min_variant_number
    ))
  }

  hap_region_variant(
    vcf = vcf,
    chrom = chrom,
    start = start,
    end = end,
    samples = samples,
    variant_type = variant_type,
    genotype_mode = genotype_mode,
    missing_genotype = missing_genotype,
    min_variant_number = min_variant_number
  )
}

build_haplotype_from_region <- function(vcf,
                                        region,
                                        samples = NULL,
                                        variant_type = NULL,
                                        genotype_mode = c("code", "string"),
                                        missing_genotype = NA_character_,
                                        min_variant_number = NULL) {
  genotype_mode <- match.arg(genotype_mode)

  vt <- retrieve_vcf(
    vcf,
    chrom = region$chrom,
    start = region$start,
    end = region$end,
    variant_type = variant_type,
    as = "VariantTrack"
  )

  sample_cols <- get_vcf_sample_columns(vt)
  stop_if_not(length(sample_cols) > 0L, "No VCF sample genotype columns were found. Re-read the VCF with `keep_genotype = TRUE`.")

  if (!is.null(samples)) {
    samples <- as.character(samples)
    missing_samples <- setdiff(samples, sample_cols)
    stop_if_not(length(missing_samples) == 0L, paste0("Samples not found in VCF: ", paste(missing_samples, collapse = ", ")))
    sample_cols <- samples
  }

  stop_if_not(nrow(vt$data) > 0L, "No variants were found in the selected region.")

  geno_long <- extract_vcf_genotype_long(
    vt$data,
    sample_cols = sample_cols,
    genotype_mode = genotype_mode,
    missing_genotype = missing_genotype
  )

  stop_if_not(nrow(geno_long) > 0L, "No genotypes were found in the selected region.")

  geno_wide <- data.table::dcast(
    geno_long,
    sample_id ~ variant_id,
    value.var = "genotype",
    fill = missing_genotype
  )

  variant_order <- vt$data[["variant_id"]]
  variant_order <- variant_order[variant_order %in% names(geno_wide)]
  stop_if_not(length(variant_order) > 0L, "No genotype columns matched the retained variants.")
  data.table::setcolorder(geno_wide, c("sample_id", variant_order))

  non_missing_counts <- geno_long[, .(
    non_missing_variant_n = sum(!genotype_missing)
  ), by = sample_id]
  geno_wide[, "non_missing_variant_n" := non_missing_counts$non_missing_variant_n[
    match(sample_id, non_missing_counts$sample_id)
  ]]

  if (is.null(min_variant_number)) {
    min_variant_number <- length(variant_order)
  }
  min_variant_number <- as.integer(min_variant_number)[1L]
  stop_if_not(!is.na(min_variant_number) && min_variant_number >= 0L, "`min_variant_number` must be a non-negative integer or NULL.")

  geno_wide <- geno_wide[non_missing_variant_n >= min_variant_number]
  stop_if_not(nrow(geno_wide) > 0L, "No samples remained after `min_variant_number` filtering.")

  geno_wide[, "hap_pattern" := do.call(paste, c(lapply(.SD, hap_value_to_string), sep = "|")), .SDcols = variant_order]

  hap_map <- geno_wide[, .(sample_n = .N, samples = paste(sample_id, collapse = ";")), by = hap_pattern]
  data.table::setorder(hap_map, -sample_n, hap_pattern)
  hap_map[, "hap_id" := paste0("Hap", seq_len(.N))]

  geno_wide <- merge(geno_wide, hap_map[, .(hap_pattern, hap_id)], by = "hap_pattern", all.x = TRUE)
  data.table::setcolorder(geno_wide, c("sample_id", "hap_id", "hap_pattern", "non_missing_variant_n", variant_order))

  hap_alleles <- geno_wide[, c("hap_id", variant_order), with = FALSE]
  hap_alleles <- unique(hap_alleles, by = "hap_id")
  hap_alleles <- merge(hap_map[, .(hap_id, sample_n, samples)], hap_alleles, by = "hap_id", all.x = TRUE)
  data.table::setorder(hap_alleles, -sample_n, hap_id)

  out <- list(
    region = region,
    variants = vt$data[],
    genotype_long = geno_long[],
    genotype_wide = geno_wide[],
    haplotypes = hap_alleles[],
    sample_haplotypes = geno_wide[, .(sample_id, hap_id, hap_pattern, non_missing_variant_n)],
    meta = list(
      genotype_mode = genotype_mode,
      missing_genotype = missing_genotype,
      min_variant_number = min_variant_number,
      sample_n = length(unique(geno_wide$sample_id)),
      variant_n = nrow(vt$data),
      haplotype_n = nrow(hap_alleles)
    )
  )

  structure(out, class = "HapVariant")
}

#' @export
print.HapVariant <- function(x, ...) {
  cat("<HapVariant>\n")
  cat("  region    : ", x$region$chrom, ":", x$region$start, "-", x$region$end, "\n", sep = "")
  cat("  variants  : ", format(nrow(x$variants), big.mark = ","), "\n", sep = "")
  cat("  samples   : ", format(x$meta$sample_n %||% 0L, big.mark = ","), "\n", sep = "")
  cat("  haplotypes: ", format(nrow(x$haplotypes), big.mark = ","), "\n", sep = "")
  invisible(x)
}

resolve_haplotype_gene_region <- function(annotation,
                                          gene_id = NULL,
                                          transcript_id = NULL,
                                          upstream = 0L,
                                          downstream = 0L,
                                          strand_aware = TRUE) {
  stop_if_not(!is.null(annotation), "`annotation` is required for gene/transcript haplotype queries.")
  locators <- sum(c(!is.null(gene_id), !is.null(transcript_id)))
  stop_if_not(locators == 1L, "Specify exactly one of `gene_id` or `transcript_id`.")

  upstream <- as.integer(upstream)[1L]
  downstream <- as.integer(downstream)[1L]
  if (is.na(upstream)) upstream <- 0L
  if (is.na(downstream)) downstream <- 0L
  stop_if_not(upstream >= 0L, "`upstream` must be a non-negative integer.")
  stop_if_not(downstream >= 0L, "`downstream` must be a non-negative integer.")

  if (!is.null(gene_id)) {
    target_id <- as.character(gene_id)[1L]
    x <- retrieve_feature(annotation, gene_id = target_id, as = "Feature")
    stop_if_not(nrow(x$genes) > 0L, "Gene ID was not found in annotation.")
    dt <- data.table::as.data.table(x$genes)
    chrom <- as.character(dt[["chrom"]][1L])
    strand <- if ("strand" %in% names(dt)) as.character(dt[["strand"]][1L]) else NA_character_
    core_start <- as.integer(min(dt[["gene_start"]], na.rm = TRUE))
    core_end <- as.integer(max(dt[["gene_end"]], na.rm = TRUE))
    locator <- "gene"
  } else {
    target_id <- as.character(transcript_id)[1L]
    x <- retrieve_feature(annotation, transcript_id = target_id, as = "Feature")
    stop_if_not(nrow(x$transcripts) > 0L, "Transcript ID was not found in annotation.")
    dt <- data.table::as.data.table(x$transcripts)
    chrom <- as.character(dt[["chrom"]][1L])
    strand <- if ("strand" %in% names(dt)) as.character(dt[["strand"]][1L]) else NA_character_
    core_start <- as.integer(min(dt[["tx_start"]], na.rm = TRUE))
    core_end <- as.integer(max(dt[["tx_end"]], na.rm = TRUE))
    locator <- "transcript"
  }

  if (isTRUE(strand_aware) && identical(strand, "-")) {
    region_start <- core_start - downstream
    region_end <- core_end + upstream
  } else {
    region_start <- core_start - upstream
    region_end <- core_end + downstream
  }

  region_start <- max(1L, as.integer(region_start))
  region_end <- max(region_start, as.integer(region_end))

  list(
    locator = locator,
    id = target_id,
    chrom = chrom,
    start = region_start,
    end = region_end,
    core_start = core_start,
    core_end = core_end,
    upstream = upstream,
    downstream = downstream,
    strand = strand,
    strand_aware = isTRUE(strand_aware)
  )
}

get_vcf_sample_columns <- function(vt) {
  meta_samples <- vt$meta$sample_names %||% character()
  meta_samples <- as.character(meta_samples)
  meta_samples <- meta_samples[meta_samples %in% names(vt$data)]
  if (length(meta_samples) > 0L) return(meta_samples)

  standard_cols <- c("chrom", "pos", "start", "end", "variant_id", "ref", "alt", "qual", "filter", "info", "variant_type", "FORMAT")
  setdiff(names(vt$data), standard_cols)
}

extract_vcf_genotype_long <- function(dt, sample_cols, genotype_mode = "code", missing_genotype = NA_character_) {
  stop_if_not(length(sample_cols) > 0L, "No sample columns were supplied.")
  id_cols <- intersect(c("chrom", "pos", "variant_id", "ref", "alt", "variant_type", "FORMAT"), names(dt))
  x <- data.table::melt(
    dt,
    id.vars = id_cols,
    measure.vars = sample_cols,
    variable.name = "sample_id",
    value.name = "genotype_raw",
    variable.factor = FALSE
  )
  if (!"FORMAT" %in% names(x)) {
    x[, "FORMAT" := "GT"]
  }
  x[, "genotype" := normalize_vcf_gt(
    genotype_raw,
    FORMAT,
    ref,
    alt,
    mode = genotype_mode,
    missing_genotype = NA_character_
  )]
  x[, "genotype_missing" := is.na(genotype)]
  x[, "genotype" := replace_missing_label(genotype, missing_genotype)]
  x[, .(
    sample_id = as.character(sample_id),
    chrom,
    pos,
    variant_id,
    ref,
    alt,
    variant_type,
    genotype,
    genotype_missing
  )]
}

normalize_vcf_gt <- function(genotype_raw, format, ref, alt, mode = c("code", "string"), missing_genotype = NA_character_) {
  mode <- match.arg(mode)
  gt <- extract_gt_field(genotype_raw, format)
  gt <- normalize_gt_separator(gt)
  gt[is_missing_gt(gt)] <- NA_character_

  if (mode == "code") {
    return(convert_gt_to_code(gt, missing_genotype = missing_genotype))
  }

  convert_gt_to_string(gt, ref, alt, missing_genotype = missing_genotype)
}

extract_gt_field <- function(x, format) {
  x <- as.character(x)
  format <- as.character(format)
  gt <- x
  has_format <- !is.na(format) & nzchar(format)
  if (any(has_format)) {
    fmt_list <- strsplit(format[has_format], ":", fixed = TRUE)
    gt_index <- vapply(fmt_list, function(z) match("GT", z), integer(1L))
    raw_list <- strsplit(x[has_format], ":", fixed = TRUE)
    gt[has_format] <- vapply(seq_along(raw_list), function(i) {
      idx <- gt_index[i]
      if (is.na(idx) || idx < 1L || idx > length(raw_list[[i]])) return(NA_character_)
      raw_list[[i]][idx]
    }, character(1L))
  }
  gt
}

normalize_gt_separator <- function(gt) {
  gt <- as.character(gt)
  gt <- gsub("/", "|", gt, fixed = TRUE)
  gt
}

is_missing_gt <- function(gt) {
  is.na(gt) | gt == "." | gt == "" | grepl("\\.", gt)
}

replace_missing_label <- function(x, missing_genotype = NA_character_) {
  x[is.na(x)] <- missing_genotype
  x
}

hap_value_to_string <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x
}

convert_gt_to_code <- function(gt, missing_genotype = NA_character_) {
  out <- rep(NA_character_, length(gt))

  for (i in seq_along(gt)) {
    g <- gt[i]
    if (is.na(g)) {
      out[i] <- missing_genotype
      next
    }

    idx <- suppressWarnings(as.integer(strsplit(g, "\\|", fixed = FALSE)[[1L]]))
    if (length(idx) == 0L || any(is.na(idx))) {
      out[i] <- missing_genotype
      next
    }

    out[i] <- if (any(idx > 0L)) "1" else "0"
  }

  replace_missing_label(out, missing_genotype)
}

convert_gt_to_string <- function(gt, ref, alt, missing_genotype = NA_character_) {
  out <- rep(NA_character_, length(gt))

  for (i in seq_along(gt)) {
    g <- gt[i]
    if (is.na(g)) {
      out[i] <- missing_genotype
      next
    }

    alt_values <- strsplit(as.character(alt[i]), ",", fixed = TRUE)[[1L]]
    alt_values <- alt_values[!is.na(alt_values) & nzchar(alt_values)]
    allele_values <- c(as.character(ref[i]), alt_values)
    idx_raw <- suppressWarnings(as.integer(strsplit(g, "\\|", fixed = FALSE)[[1L]]))

    if (length(idx_raw) == 0L || any(is.na(idx_raw)) || any(idx_raw < 0L)) {
      out[i] <- missing_genotype
      next
    }

    # Display one compact allele label per variant. If any ALT allele is present,
    # show the first observed ALT allele; otherwise show the REF allele.
    display_index <- if (any(idx_raw > 0L)) idx_raw[idx_raw > 0L][1L] + 1L else 1L
    if (display_index < 1L || display_index > length(allele_values)) {
      out[i] <- missing_genotype
      next
    }

    value <- format_hap_allele(allele_values[display_index])
    if (is.na(value) || !nzchar(value)) {
      out[i] <- missing_genotype
    } else {
      out[i] <- value
    }
  }

  replace_missing_label(out, missing_genotype)
}

format_hap_allele <- function(x) {
  x <- as.character(x)
  out <- x
  invalid <- is.na(out) | out == "" | out == "."
  is_compact_indel <- !invalid & grepl("^i[0-9]+$", out)
  allele_len <- nchar(out, type = "chars", allowNA = TRUE, keepNA = TRUE)
  is_indel_like <- !invalid & !is_compact_indel & !is.na(allele_len) & allele_len != 1L
  out[is_indel_like] <- paste0("i", allele_len[is_indel_like])
  out[invalid] <- NA_character_
  out
}

format_hap_ref_alt <- function(ref, alt) {
  ref <- as.character(ref)
  alt <- as.character(alt)
  n <- max(length(ref), length(alt))
  ref <- rep(ref, length.out = n)
  alt <- rep(alt, length.out = n)

  vapply(seq_len(n), function(i) {
    ref_i <- format_hap_allele(ref[i])
    alt_i <- strsplit(alt[i], ",", fixed = TRUE)[[1L]]
    alt_i <- alt_i[!is.na(alt_i) & nzchar(alt_i)]
    alt_i <- format_hap_allele(alt_i)
    alt_i <- alt_i[!is.na(alt_i) & nzchar(alt_i)]
    if (is.na(ref_i) || !nzchar(ref_i)) ref_i <- "NA"
    if (length(alt_i) == 0L) alt_i <- "NA"
    paste0(ref_i, "/", paste(alt_i, collapse = ","))
  }, character(1L))
}
