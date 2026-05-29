# Author: Rensc
# Date: 2026-05-29
# Version: 0.1.0
# Function: Build haplotype tables from VCF variants
# Input: VariantTrack objects, VCF files, and genomic locators
# Output: HapVariant objects

#' Build haplotypes from variants in a gene or genomic region
#'
#' @description
#' Extracts variants from a VCF/VariantTrack object and converts sample genotype
#' profiles into haplotypes. The function supports region-based queries and
#' annotation-guided gene/transcript queries.
#'
#' @param vcf A VariantTrack object or VCF file path.
#' @param annotation Optional gene annotation object used for `gene_id` or `transcript_id` queries.
#' @param gene_id Optional gene ID.
#' @param transcript_id Optional transcript ID.
#' @param chrom Optional chromosome name.
#' @param start Optional region start in 1-based closed coordinates.
#' @param end Optional region end in 1-based closed coordinates.
#' @param samples Optional sample names to keep.
#' @param variant_type Optional variant types to keep.
#' @param genotype_mode Genotype representation: `code` keeps GT codes, `allele` converts GT to REF/ALT alleles.
#' @param missing_genotype Missing genotype symbol.
#' @param min_variant_non_missing Minimum number of non-missing samples required for a variant.
#' @param min_hap_samples Minimum number of samples required to keep a haplotype group.
#' @return A HapVariant object.
#' @export
hap_variant <- function(vcf,
                        annotation = NULL,
                        gene_id = NULL,
                        transcript_id = NULL,
                        chrom = NULL,
                        start = NULL,
                        end = NULL,
                        samples = NULL,
                        variant_type = NULL,
                        genotype_mode = c("code", "allele"),
                        missing_genotype = "./.",
                        min_variant_non_missing = 1L,
                        min_hap_samples = 1L) {
  genotype_mode <- match.arg(genotype_mode)

  region <- resolve_haplotype_region(
    annotation = annotation,
    gene_id = gene_id,
    transcript_id = transcript_id,
    chrom = chrom,
    start = start,
    end = end
  )

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

  geno_long <- extract_vcf_genotype_long(
    vt$data,
    sample_cols = sample_cols,
    genotype_mode = genotype_mode,
    missing_genotype = missing_genotype
  )

  if (nrow(geno_long) > 0L) {
    non_missing <- geno_long[, .(non_missing_n = sum(!is.na(genotype) & genotype != missing_genotype)), by = variant_id]
    keep_vars <- non_missing[non_missing_n >= as.integer(min_variant_non_missing), variant_id]
    geno_long <- geno_long[variant_id %in% keep_vars]
    vt$data <- vt$data[variant_id %in% keep_vars]
  }

  stop_if_not(nrow(vt$data) > 0L, "No variants remained after filtering.")

  geno_wide <- data.table::dcast(
    geno_long,
    sample_id ~ variant_id,
    value.var = "genotype",
    fill = missing_genotype
  )

  variant_order <- vt$data[["variant_id"]]
  variant_order <- variant_order[variant_order %in% names(geno_wide)]
  data.table::setcolorder(geno_wide, c("sample_id", variant_order))

  geno_wide[, "hap_pattern" := do.call(paste, c(.SD, sep = "|")), .SDcols = variant_order]
  hap_map <- geno_wide[, .(sample_n = .N, samples = paste(sample_id, collapse = ";")), by = hap_pattern]
  data.table::setorder(hap_map, -sample_n, hap_pattern)
  hap_map[, "hap_id" := paste0("Hap", seq_len(.N))]
  hap_map <- hap_map[sample_n >= as.integer(min_hap_samples)]

  geno_wide <- merge(geno_wide, hap_map[, .(hap_pattern, hap_id)], by = "hap_pattern", all.x = FALSE)
  data.table::setcolorder(geno_wide, c("sample_id", "hap_id", "hap_pattern", variant_order))

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
    sample_haplotypes = geno_wide[, .(sample_id, hap_id, hap_pattern)],
    meta = list(
      genotype_mode = genotype_mode,
      missing_genotype = missing_genotype,
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

resolve_haplotype_region <- function(annotation = NULL,
                                     gene_id = NULL,
                                     transcript_id = NULL,
                                     chrom = NULL,
                                     start = NULL,
                                     end = NULL) {
  locators <- sum(c(!is.null(gene_id), !is.null(transcript_id), !is.null(chrom) || !is.null(start) || !is.null(end)))
  stop_if_not(locators == 1L, "Specify exactly one locator: `gene_id`, `transcript_id`, or `chrom` + `start` + `end`.")

  if (!is.null(gene_id)) {
    stop_if_not(!is.null(annotation), "`annotation` is required when `gene_id` is used.")
    x <- retrieve_feature(annotation, gene_id = as.character(gene_id)[1L], as = "Feature")
    stop_if_not(nrow(x$genes) > 0L, "Gene ID was not found in annotation.")
    return(list(
      locator = "gene",
      id = as.character(gene_id)[1L],
      chrom = as.character(x$genes$chrom[1L]),
      start = as.integer(min(x$genes$gene_start, na.rm = TRUE)),
      end = as.integer(max(x$genes$gene_end, na.rm = TRUE))
    ))
  }

  if (!is.null(transcript_id)) {
    stop_if_not(!is.null(annotation), "`annotation` is required when `transcript_id` is used.")
    x <- retrieve_feature(annotation, transcript_id = as.character(transcript_id)[1L], as = "Feature")
    stop_if_not(nrow(x$transcripts) > 0L, "Transcript ID was not found in annotation.")
    return(list(
      locator = "transcript",
      id = as.character(transcript_id)[1L],
      chrom = as.character(x$transcripts$chrom[1L]),
      start = as.integer(min(x$transcripts$tx_start, na.rm = TRUE)),
      end = as.integer(max(x$transcripts$tx_end, na.rm = TRUE))
    ))
  }

  check_region(chrom, start, end)
  list(
    locator = "region",
    id = paste0(as.character(chrom)[1L], ":", as.integer(start)[1L], "-", as.integer(end)[1L]),
    chrom = as.character(chrom)[1L],
    start = as.integer(start)[1L],
    end = as.integer(end)[1L]
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

extract_vcf_genotype_long <- function(dt, sample_cols, genotype_mode = "code", missing_genotype = "./.") {
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
  x[, "genotype" := normalize_vcf_gt(genotype_raw, FORMAT, ref, alt, mode = genotype_mode, missing_genotype = missing_genotype)]
  x[, .(sample_id = as.character(sample_id), chrom, pos, variant_id, ref, alt, variant_type, genotype)]
}

normalize_vcf_gt <- function(genotype_raw, format, ref, alt, mode = "code", missing_genotype = "./.") {
  gt <- extract_gt_field(genotype_raw, format)
  gt[is.na(gt) | gt == "." | gt == "" | grepl("\\.", gt)] <- missing_genotype
  if (mode == "code") return(gt)
  convert_gt_to_alleles(gt, ref, alt, missing_genotype = missing_genotype)
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

convert_gt_to_alleles <- function(gt, ref, alt, missing_genotype = "./.") {
  out <- gt
  for (i in seq_along(gt)) {
    g <- gt[i]
    if (is.na(g) || g == missing_genotype) {
      out[i] <- missing_genotype
      next
    }
    sep <- if (grepl("\\|", g)) "|" else "/"
    alleles <- c(as.character(ref[i]), strsplit(as.character(alt[i]), ",", fixed = TRUE)[[1L]])
    idx <- suppressWarnings(as.integer(strsplit(g, "[/|]")[[1L]])) + 1L
    if (any(is.na(idx)) || any(idx < 1L) || any(idx > length(alleles))) {
      out[i] <- missing_genotype
    } else {
      out[i] <- paste(alleles[idx], collapse = sep)
    }
  }
  out
}
