# Author: Rensc
# Date: 2026-11-09
# Version: dev001
# Function: Recover raw genotype states and identify exact redundancy blocks
# Input: HapVariant objects and optional sample-by-variant ALT dosage matrices
# Output: Aligned genotypes, raw pattern keys, and reversible variant mappings

refine_genotype_data <- function(hap, genotypes = NULL, protect_variants = character()) {
  stop_if_not(inherits(hap, "HapVariant"), "`hap` must be a HapVariant object.")
  sites <- c("chrom", "pos", "ref", "alt")
  stop_if_not(
    is.data.frame(hap$variants) && nrow(hap$variants) > 0L &&
      !anyDuplicated(names(hap$variants)) && all(c(sites, "variant_id") %in% names(hap$variants)),
    "`hap$variants` must contain unique columns for variant IDs and genomic sites."
  )
  variants <- data.table::as.data.table(hap$variants)
  stop_if_not(refine_gt_valid_ids(variants[["variant_id"]]), "Variant IDs must be non-empty and unique.")
  pos <- variants[["pos"]]
  stop_if_not(
    is.numeric(pos) && all(is.finite(pos) & pos > 0 & pos == floor(pos)) &&
      all(vapply(c("chrom", "ref", "alt"), function(field) {
        x <- variants[[field]]
        is.character(x) && !anyNA(x) && all(nzchar(x))
      }, logical(1L))),
    "Genomic sites require chromosome names, positive integer positions, REF and ALT alleles."
  )
  stop_if_not(!anyDuplicated(variants, by = sites), "Genomic site keys (chrom, pos, ref, alt) must be unique.")
  data.table::setorderv(variants, c(sites, "variant_id"))
  variant_ids <- variants[["variant_id"]]
  stop_if_not(
    is.character(protect_variants) && all(protect_variants %in% variant_ids),
    "`protect_variants` must contain only known variant IDs."
  )

  samples <- hap$sample_haplotypes
  stop_if_not(
    is.data.frame(samples) && refine_gt_valid_ids(samples[["sample_id"]]),
    "`hap$sample_haplotypes` must contain non-empty, unique sample IDs."
  )
  hap_ids <- samples[["hap_id"]]
  stop_if_not(
    is.character(hap_ids) && !anyNA(hap_ids) && length(hap_ids) == nrow(samples) && all(nzchar(hap_ids)),
    "Each sample must have an original haplotype ID."
  )
  sample_order <- order(samples[["sample_id"]], method = "radix")
  sample_ids <- samples[["sample_id"]][sample_order]
  original_hap_id <- stats::setNames(hap_ids[sample_order], sample_ids)
  matrix_names <- list(sample_ids, variant_ids)

  ref <- toupper(variants[["ref"]])
  alt <- toupper(variants[["alt"]])
  biallelic <- grepl("^[ACGTN]+$", ref) & grepl("^[ACGTN]+$", alt) & ref != alt
  unsupported <- matrix(FALSE, length(sample_ids), length(variant_ids), dimnames = matrix_names)
  pattern <- stats::setNames(character(length(sample_ids)), sample_ids)
  if (is.null(genotypes)) {
    source <- "vcf_gt"
    stop_if_not(
      all(sample_ids %in% names(variants)),
      "The raw GT fields are unavailable for some samples. Supply a named `genotypes` dosage matrix."
    )
    format <- variants[["FORMAT"]] %||% rep("GT", nrow(variants))
    has_gt <- !is.na(format) & grepl("(^|:)GT(:|$)", format)
    stop_if_not(any(has_gt), "No raw GT fields were found in FORMAT. Supply a named `genotypes` dosage matrix.")
    format[!has_gt] <- "."
    dosage <- matrix(NA_real_, length(sample_ids), length(variant_ids), dimnames = matrix_names)
    for (i in seq_along(sample_ids)) {
      gt <- extract_gt_field(variants[[sample_ids[i]]], format)
      missing <- is.na(gt) | gt == "" | grepl("^\\.([/|]\\.)*$", gt) |
        (grepl("^[01.][/|][01.]$", gt) & grepl("\\.", gt))
      valid <- biallelic & grepl("^[01][/|][01]$", gt)
      unsupported[i, ] <- !missing & !valid
      dosage[i, valid] <- as.numeric(substr(gt[valid], 1L, 1L)) + as.numeric(substr(gt[valid], 3L, 3L))
      pattern[i] <- refine_gt_signature(dosage[i, ], unsupported[i, ], gt)
    }
  } else {
    source <- "dosage_matrix"
    stop_if_not(is.matrix(genotypes) && is.numeric(genotypes), "`genotypes` must be a numeric dosage matrix.")
    stop_if_not(
      refine_gt_valid_ids(rownames(genotypes)) && setequal(rownames(genotypes), sample_ids),
      "`genotypes` row names must match all sample IDs exactly once."
    )
    stop_if_not(
      refine_gt_valid_ids(colnames(genotypes)) && setequal(colnames(genotypes), variant_ids),
      "`genotypes` column names must match all variant IDs exactly once."
    )
    stop_if_not(
      !any(is.nan(genotypes)) && all(is.na(genotypes) | genotypes %in% 0:2),
      "`genotypes` values must be 0, 1, 2, or NA."
    )
    dosage <- genotypes[sample_ids, variant_ids, drop = FALSE]
    stop_if_not(
      all(is.na(dosage[, !biallelic, drop = FALSE])),
      "`genotypes` can supply dosages only for sequence-resolved biallelic variants."
    )
    storage.mode(dosage) <- "double"
    dimnames(dosage) <- matrix_names
    pattern[] <- vapply(seq_along(sample_ids), function(i) refine_gt_signature(dosage[i, ]), character(1L))
  }

  non_missing_n <- colSums(!is.na(dosage))
  state_n <- vapply(seq_along(variant_ids), function(j) {
    observed <- dosage[, j]
    length(unique(observed[!is.na(observed)]))
  }, integer(1L))
  unsupported_n <- colSums(unsupported)
  supported <- biallelic & unsupported_n == 0L & state_n >= 2L
  representatives <- seq_along(variant_ids)
  eligible <- which(supported)
  signatures <- vapply(eligible, function(j) {
    column <- dosage[, j]
    # Heterozygotes are invariant under allele reversal.
    orientation <- which(column %in% c(0, 2))[1L]
    if (column[orientation] == 2) column <- 2 - column
    refine_gt_signature(column)
  }, character(1L))
  representatives[eligible] <- eligible[match(signatures, signatures)]
  flipped <- logical(length(variant_ids))
  for (j in which(representatives != seq_along(variant_ids))) {
    column <- dosage[, j]
    representative <- dosage[, representatives[j]]
    flipped[j] <- identical(2 - column, representative)
    stop_if_not(flipped[j] || identical(column, representative), "Exact genotype block verification failed.")
  }

  representative_cols <- unique(representatives)
  block_index <- match(representatives, representative_cols)
  block_n <- tabulate(block_index, nbins = length(representative_cols))
  block_ids <- sprintf("rV%03d_%d", seq_along(representative_cols), block_n)
  protected <- block_index %in% block_index[variant_ids %in% protect_variants]
  variant_map <- data.table::data.table(
    variant_id = variant_ids,
    block_id = block_ids[block_index],
    representative_id = variant_ids[representatives],
    flipped = flipped,
    block_n = block_n[block_index],
    protected = protected,
    non_missing_n = non_missing_n,
    state_n = state_n,
    unsupported_n = unsupported_n
  )
  blocks <- dosage[, representative_cols, drop = FALSE]
  colnames(blocks) <- block_ids

  list(
    variants = variants, dosage = dosage, blocks = blocks, variant_map = variant_map,
    pattern = pattern, original_hap_id = original_hap_id, unsupported = unsupported,
    block_supported = stats::setNames(supported[representative_cols], block_ids),
    source = source
  )
}

refine_gt_valid_ids <- function(ids) {
  is.character(ids) && length(ids) > 0L && !anyNA(ids) && all(nzchar(ids)) && !anyDuplicated(ids)
}

refine_gt_signature <- function(dosage, unsupported = NULL, gt = NULL) {
  states <- ifelse(is.na(dosage), "M", as.character(dosage))
  if (any(unsupported)) {
    raw <- enc2utf8(gt[unsupported])
    # Length-prefix raw calls to keep signatures unique even for malformed GTs.
    states[unsupported] <- paste0("U", nchar(raw, type = "bytes"), ":", raw)
  }
  paste0(states, collapse = "")
}
