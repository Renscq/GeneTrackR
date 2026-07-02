# Author: Rensc
# Date: 2026-07-02
# Version: 0.2.0
# Function: LD block calculation from VCF genotypes
# Input: VariantTrack objects or VCF paths
# Output: LDTrack objects

#' Compute a linkage disequilibrium block from VCF genotypes
#'
#' @description
#' Calculates pairwise LD for variants in a genomic region. Genotypes are parsed
#' from the sample columns preserved by `read_vcf(..., keep_genotype = TRUE)`.
#' The current implementation uses genotype dosage, where 0 is homozygous
#' reference, 1 is heterozygous/non-reference carrier, and 2 is homozygous
#' alternate for diploid biallelic records. Multi-allelic non-reference alleles
#' are collapsed into total non-reference dosage.
#'
#' @param vcf A `VariantTrack` object returned by `read_vcf()` or a VCF path.
#' @param chrom Chromosome. If NULL, a single-chromosome in-memory VCF uses its
#' full chromosome range.
#' @param start Region start. If NULL, the minimum retained variant position is used.
#' @param end Region end. If NULL, the maximum retained variant position is used.
#' @param variant_type One of `both`, `snp`, or `ind`. `ind` keeps INS/DEL/MNV
#' records.
#' @param method LD value used as the primary `ld` column. Currently `r2` and
#' `Dprime` are supported. Dprime is a dosage-based approximation for unphased VCFs.
#' @param samples Optional sample names to keep.
#' @param min_pair_samples Minimum number of paired non-missing samples required
#' for a pairwise estimate.
#' @param ploidy Ploidy used for dosage-based allele-frequency and Dprime
#' approximation. Default is 2.
#' @param keep_genotype_matrix Logical. Whether to store the variant-by-sample
#' dosage matrix in the returned object.
#' @param verbose Logical. Whether to print progress messages.
#' @return An `LDTrack` object.
#' @export
compute_ld_block <- function(vcf,
                             chrom = NULL,
                             start = NULL,
                             end = NULL,
                             variant_type = c("both", "snp", "ind"),
                             method = c("r2", "Dprime"),
                             samples = NULL,
                             min_pair_samples = 3L,
                             ploidy = 2L,
                             keep_genotype_matrix = FALSE,
                             verbose = TRUE) {
  variant_type <- match.arg(variant_type)
  method <- match.arg(method)
  verbose <- isTRUE(verbose)
  min_pair_samples <- as.integer(min_pair_samples)[1L]
  stop_if_not(!is.na(min_pair_samples) && min_pair_samples >= 2L,
              "`min_pair_samples` must be an integer >= 2.")
  ploidy <- as.numeric(ploidy)[1L]
  stop_if_not(is.finite(ploidy) && ploidy > 0,
              "`ploidy` must be a positive numeric value.")

  region <- resolve_ld_region(vcf, chrom = chrom, start = start, end = end)

  vt <- retrieve_vcf(
    vcf,
    chrom = region$chrom,
    start = region$start,
    end = region$end,
    keep_genotype = TRUE,
    as = "VariantTrack",
    verbose = FALSE
  )

  vt$data <- filter_ld_variant_type(vt$data, variant_type = variant_type)
  stop_if_not(nrow(vt$data) >= 2L, "At least two variants are required for LD calculation.")

  sample_cols <- get_vcf_sample_columns(vt)
  stop_if_not(length(sample_cols) > 0L,
              "No VCF sample genotype columns were found. Re-read the VCF with `keep_genotype = TRUE`.")

  if (!is.null(samples)) {
    samples <- as.character(samples)
    missing_samples <- setdiff(samples, sample_cols)
    stop_if_not(length(missing_samples) == 0L,
                paste0("Samples not found in VCF: ", paste(missing_samples, collapse = ", ")))
    sample_cols <- samples
  }

  data.table::setorderv(vt$data, c("chrom", "pos", "variant_id"))
  vt$data[, "variant_index" := seq_len(.N)]
  variants <- vt$data[, intersect(
    c("variant_index", "chrom", "pos", "variant_id", "ref", "alt", "variant_type"),
    names(vt$data)
  ), with = FALSE]

  if (verbose && nrow(variants) > 2000L) {
    warning(
      "The selected region contains more than 2,000 variants; pairwise LD output may be large.",
      call. = FALSE
    )
  }

  if (verbose) {
    message("[GeneTrackR] Computing LD for ", nrow(variants), " variants and ", length(sample_cols), " samples.")
  }

  dosage <- extract_ld_dosage_matrix(vt$data, sample_cols = sample_cols)
  pair_dt <- compute_ld_pair_table(
    dosage = dosage,
    variants = variants,
    method = method,
    min_pair_samples = min_pair_samples,
    ploidy = ploidy
  )
  ld_matrix <- build_ld_matrix(pair_dt, variants = variants)

  meta <- list(
    method = method,
    variant_type = variant_type,
    sample_n = length(sample_cols),
    min_pair_samples = min_pair_samples,
    ploidy = ploidy,
    dosage_definition = "number of non-reference allele copies parsed from GT",
    dprime_note = "Dprime is approximated from genotype dosage for unphased VCFs. Use phased haplotypes for exact haplotype-based Dprime."
  )

  LDTrack(
    data = pair_dt,
    matrix = ld_matrix,
    variants = variants,
    region = region,
    genotype = if (isTRUE(keep_genotype_matrix)) dosage else NULL,
    meta = meta
  )
}

resolve_ld_region <- function(vcf, chrom = NULL, start = NULL, end = NULL) {
  has_any_region <- !is.null(chrom) || !is.null(start) || !is.null(end)
  if (!is.null(start) || !is.null(end)) {
    stop_if_not(!is.null(chrom), "`chrom` is required when `start` or `end` is supplied.")
    stop_if_not(!is.null(start) && !is.null(end), "Both `start` and `end` are required for LD region queries.")
  }

  if (inherits(vcf, "VariantTrack") && is_lazy_variant_track(vcf)) {
    stop_if_not(!is.null(chrom) && !is.null(start) && !is.null(end),
                "Lazy VariantTrack LD queries require `chrom`, `start`, and `end`.")
    return(list(
      chrom = as.character(chrom)[1L],
      start = as.integer(start)[1L],
      end = as.integer(end)[1L]
    ))
  }

  if (is.character(vcf) && length(vcf) == 1L && file.exists(vcf)) {
    if (isTRUE(has_any_region)) {
      stop_if_not(!is.null(chrom) && !is.null(start) && !is.null(end),
                  "VCF path LD queries require complete `chrom`, `start`, and `end` unless the whole file is read.")
      return(list(
        chrom = as.character(chrom)[1L],
        start = as.integer(start)[1L],
        end = as.integer(end)[1L]
      ))
    }
    vt <- read_vcf(vcf, keep_genotype = TRUE, mode = "memory", verbose = FALSE)
  } else {
    stop_if_not(inherits(vcf, "VariantTrack"), "`vcf` must be a VariantTrack object or a VCF file path.")
    vt <- vcf
  }

  dt <- data.table::as.data.table(vt$data)
  stop_if_not(nrow(dt) > 0L, "Cannot infer LD region from an empty VCF. Please supply `chrom`, `start`, and `end`.")

  if (!is.null(chrom)) {
    chrom_value <- as.character(chrom)[1L]
    dt <- dt[as.character(dt[["chrom"]]) == chrom_value]
    stop_if_not(nrow(dt) > 0L, "No variants were found on the requested chromosome.")
  } else {
    chrom_values <- unique(as.character(dt[["chrom"]]))
    chrom_values <- chrom_values[!is.na(chrom_values) & nzchar(chrom_values)]
    stop_if_not(length(chrom_values) == 1L,
                "VCF contains multiple chromosomes. Please supply `chrom`, `start`, and `end` for LD calculation.")
    chrom_value <- chrom_values[1L]
  }

  start_value <- if (is.null(start)) min(as.integer(dt[["pos"]]), na.rm = TRUE) else as.integer(start)[1L]
  end_value <- if (is.null(end)) max(as.integer(dt[["pos"]]), na.rm = TRUE) else as.integer(end)[1L]
  stop_if_not(!is.na(start_value) && !is.na(end_value) && start_value <= end_value,
              "Invalid LD region: `start` must be <= `end`.")

  list(chrom = chrom_value, start = start_value, end = end_value)
}

filter_ld_variant_type <- function(dt, variant_type = c("both", "snp", "ind")) {
  variant_type <- match.arg(variant_type)
  x <- data.table::copy(data.table::as.data.table(dt))
  if (variant_type == "both" || nrow(x) == 0L) return(x[])
  type <- toupper(as.character(x[["variant_type"]]))
  if (variant_type == "snp") {
    return(x[type == "SNP"][])
  }
  x[type %in% c("INS", "DEL", "MNV", "INDEL")][]
}

extract_ld_dosage_matrix <- function(dt, sample_cols) {
  stop_if_not(length(sample_cols) > 0L, "No sample genotype columns were supplied.")
  id_cols <- intersect(c("variant_index", "variant_id", "FORMAT"), names(dt))
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
  x[, "dosage" := vcf_gt_to_alt_dosage(genotype_raw, FORMAT)]
  wide <- data.table::dcast(
    x,
    variant_index + variant_id ~ sample_id,
    value.var = "dosage"
  )
  data.table::setorderv(wide, "variant_index")
  sample_cols <- sample_cols[sample_cols %in% names(wide)]
  mat <- as.matrix(wide[, sample_cols, with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- as.character(wide[["variant_id"]])
  mat
}

vcf_gt_to_alt_dosage <- function(genotype_raw, format) {
  gt <- extract_gt_field(genotype_raw, format)
  gt <- normalize_gt_separator(gt)
  gt[is_missing_gt(gt)] <- NA_character_
  vapply(gt, function(z) {
    if (is.na(z) || !nzchar(z)) return(NA_real_)
    alleles <- strsplit(z, "|", fixed = TRUE)[[1L]]
    alleles <- alleles[!is.na(alleles) & alleles != "." & alleles != ""]
    if (length(alleles) == 0L) return(NA_real_)
    sum(alleles != "0")
  }, numeric(1L), USE.NAMES = FALSE)
}

compute_ld_pair_table <- function(dosage, variants, method, min_pair_samples = 3L, ploidy = 2L) {
  n_var <- nrow(variants)
  if (n_var < 2L) return(data.table::data.table())

  out <- vector("list", n_var - 1L)
  for (i in seq_len(n_var - 1L)) {
    js <- (i + 1L):n_var
    xi <- dosage[i, ]
    block <- lapply(js, function(j) {
      stats <- compute_ld_pair_stats(
        x = xi,
        y = dosage[j, ],
        min_pair_samples = min_pair_samples,
        ploidy = ploidy
      )
      data.table::data.table(
        variant_i = as.character(variants[["variant_id"]][i]),
        variant_j = as.character(variants[["variant_id"]][j]),
        index_i = as.integer(variants[["variant_index"]][i]),
        index_j = as.integer(variants[["variant_index"]][j]),
        pos_i = as.integer(variants[["pos"]][i]),
        pos_j = as.integer(variants[["pos"]][j]),
        distance_bp = as.integer(abs(as.integer(variants[["pos"]][j]) - as.integer(variants[["pos"]][i]))),
        n_samples = as.integer(stats$n_samples),
        r = as.numeric(stats$r),
        r2 = as.numeric(stats$r2),
        D = as.numeric(stats$D),
        Dprime_signed = as.numeric(stats$Dprime_signed),
        Dprime = as.numeric(stats$Dprime),
        p_i = as.numeric(stats$p_i),
        p_j = as.numeric(stats$p_j)
      )
    })
    out[[i]] <- data.table::rbindlist(block, fill = TRUE)
  }

  pair_dt <- data.table::rbindlist(out, fill = TRUE)
  pair_dt[, "ld" := if (identical(method, "r2")) as.numeric(r2) else as.numeric(Dprime)]
  pair_dt[, "method" := method]
  pair_dt[]
}

compute_ld_pair_stats <- function(x, y, min_pair_samples = 3L, ploidy = 2L) {
  ok <- !is.na(x) & !is.na(y)
  n <- sum(ok)
  empty <- list(
    n_samples = n,
    r = NA_real_,
    r2 = NA_real_,
    D = NA_real_,
    Dprime_signed = NA_real_,
    Dprime = NA_real_,
    p_i = NA_real_,
    p_j = NA_real_
  )
  if (n < min_pair_samples) return(empty)

  x <- as.numeric(x[ok])
  y <- as.numeric(y[ok])
  if (stats::var(x) <= 0 || stats::var(y) <= 0) return(empty)

  r <- suppressWarnings(stats::cor(x, y))
  r2 <- if (is.na(r)) NA_real_ else max(0, min(1, r * r))

  p_i <- mean(x, na.rm = TRUE) / ploidy
  p_j <- mean(y, na.rm = TRUE) / ploidy
  p_i <- max(0, min(1, p_i))
  p_j <- max(0, min(1, p_j))

  d_cov <- mean((x - mean(x)) * (y - mean(y)))
  D <- d_cov / ploidy
  Dmax <- if (is.na(D)) {
    NA_real_
  } else if (D >= 0) {
    min(p_i * (1 - p_j), (1 - p_i) * p_j)
  } else {
    min(p_i * p_j, (1 - p_i) * (1 - p_j))
  }
  Dprime_signed <- if (is.na(Dmax) || Dmax <= 0) NA_real_ else D / Dmax
  if (!is.na(Dprime_signed)) {
    Dprime_signed <- max(-1, min(1, Dprime_signed))
  }
  Dprime <- abs(Dprime_signed)

  list(
    n_samples = n,
    r = r,
    r2 = r2,
    D = D,
    Dprime_signed = Dprime_signed,
    Dprime = Dprime,
    p_i = p_i,
    p_j = p_j
  )
}

build_ld_matrix <- function(pair_dt, variants) {
  n_var <- nrow(variants)
  ids <- as.character(variants[["variant_id"]])
  mat <- matrix(NA_real_, nrow = n_var, ncol = n_var, dimnames = list(ids, ids))
  diag(mat) <- 1
  if (nrow(pair_dt) == 0L) return(mat)
  ii <- as.integer(pair_dt[["index_i"]])
  jj <- as.integer(pair_dt[["index_j"]])
  vals <- as.numeric(pair_dt[["ld"]])
  mat[cbind(ii, jj)] <- vals
  mat[cbind(jj, ii)] <- vals
  mat
}

