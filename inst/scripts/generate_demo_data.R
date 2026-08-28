# Author: Rensc
# Date: 2026-08-23
# Version: dev008
# Function: Generate the deterministic GeneTrackR demo dataset from canonical model tables
# Input: inst/scripts/demo_model/*.tsv
# Output: inst/extdata/gtr_demo_* example input files

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

read_tsv <- function(file) {
  utils::read.delim(
    file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

parse_int_list <- function(x) {
  as.integer(strsplit(as.character(x), ",", fixed = TRUE)[[1L]])
}

get_exons <- function(tx_row) {
  starts <- parse_int_list(tx_row$exon_starts)
  ends <- parse_int_list(tx_row$exon_ends)
  if (length(starts) != length(ends)) {
    stop("Exon start/end counts do not match for ", tx_row$transcript_id, call. = FALSE)
  }
  data.frame(start = starts, end = ends, stringsAsFactors = FALSE)
}

build_gene_table <- function(transcripts) {
  gene_ids <- unique(transcripts$gene_id)
  out <- lapply(gene_ids, function(gene_id) {
    x <- transcripts[transcripts$gene_id == gene_id, , drop = FALSE]
    data.frame(
      gene_id = gene_id,
      chrom = x$chrom[1L],
      strand = x$strand[1L],
      gene_start = min(x$tx_start),
      gene_end = max(x$tx_end),
      gene_type = x$gene_type[1L],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out[order(out$chrom, out$gene_start, out$gene_id), , drop = FALSE]
}

build_subfeatures <- function(tx_row) {
  exons <- get_exons(tx_row)
  strand <- as.character(tx_row$strand)
  tx_order <- if (strand == "+") seq_len(nrow(exons)) else rev(seq_len(nrow(exons)))
  exon_number <- integer(nrow(exons))
  exon_number[tx_order] <- seq_along(tx_order)

  has_cds <- !is.na(tx_row$cds_start) && !is.na(tx_row$cds_end) &&
    nzchar(as.character(tx_row$cds_start)) && nzchar(as.character(tx_row$cds_end))
  cds_start <- if (has_cds) as.integer(tx_row$cds_start) else NA_integer_
  cds_end <- if (has_cds) as.integer(tx_row$cds_end) else NA_integer_

  rows <- list()
  k <- 0L
  for (i in seq_len(nrow(exons))) {
    start <- exons$start[i]
    end <- exons$end[i]
    k <- k + 1L
    rows[[k]] <- data.frame(type = "exon", start = start, end = end, phase = ".", exon_number = exon_number[i])

    if (!has_cds) next

    if (start < cds_start) {
      utr_start <- start
      utr_end <- min(end, cds_start - 1L)
      if (utr_start <= utr_end) {
        utr_type <- if (strand == "+") "five_prime_UTR" else "three_prime_UTR"
        k <- k + 1L
        rows[[k]] <- data.frame(type = utr_type, start = utr_start, end = utr_end, phase = ".", exon_number = exon_number[i])
      }
    }

    cds_seg_start <- max(start, cds_start)
    cds_seg_end <- min(end, cds_end)
    if (cds_seg_start <= cds_seg_end) {
      k <- k + 1L
      rows[[k]] <- data.frame(type = "CDS", start = cds_seg_start, end = cds_seg_end, phase = NA_character_, exon_number = exon_number[i])
    }

    if (end > cds_end) {
      utr_start <- max(start, cds_end + 1L)
      utr_end <- end
      if (utr_start <= utr_end) {
        utr_type <- if (strand == "+") "three_prime_UTR" else "five_prime_UTR"
        k <- k + 1L
        rows[[k]] <- data.frame(type = utr_type, start = utr_start, end = utr_end, phase = ".", exon_number = exon_number[i])
      }
    }
  }

  out <- do.call(rbind, rows)
  cds_idx <- which(out$type == "CDS")
  if (length(cds_idx) > 0L) {
    cds_order <- cds_idx[order(out$start[cds_idx], decreasing = strand == "-")]
    cumulative <- 0L
    for (idx in cds_order) {
      out$phase[idx] <- as.character((3L - (cumulative %% 3L)) %% 3L)
      cumulative <- cumulative + out$end[idx] - out$start[idx] + 1L
    }
  }
  out[order(out$start, match(out$type, c("exon", "five_prime_UTR", "CDS", "three_prime_UTR")), out$end), , drop = FALSE]
}

write_genepred_ext <- function(transcripts, file) {
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  for (i in seq_len(nrow(transcripts))) {
    x <- transcripts[i, , drop = FALSE]
    exons <- get_exons(x)
    tx_start <- as.integer(x$tx_start) - 1L
    tx_end <- as.integer(x$tx_end)
    has_cds <- !is.na(x$cds_start) && nzchar(as.character(x$cds_start))
    cds_start <- if (has_cds) as.integer(x$cds_start) - 1L else tx_start
    cds_end <- if (has_cds) as.integer(x$cds_end) else tx_start
    cds_stat <- if (has_cds) "complete" else "none"
    exon_starts <- paste0(paste(exons$start - 1L, collapse = ","), ",")
    exon_ends <- paste0(paste(exons$end, collapse = ","), ",")
    exon_frames <- paste0(paste(rep("-1", nrow(exons)), collapse = ","), ",")
    fields <- c(
      x$transcript_id, x$chrom, x$strand, tx_start, tx_end, cds_start, cds_end,
      nrow(exons), exon_starts, exon_ends, 0, x$gene_id, cds_stat, cds_stat, exon_frames
    )
    writeLines(paste(fields, collapse = "\t"), con)
  }
}

write_gff3 <- function(genes, transcripts, file) {
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c("##gff-version 3", "##sequence-region chr1 1 25000000", "##sequence-region chr2 1 18000000"), con)
  for (i in seq_len(nrow(genes))) {
    gene <- genes[i, , drop = FALSE]
    gene_attr <- paste0("ID=", gene$gene_id, ";Name=", gene$gene_id, ";gene_biotype=", gene$gene_type)
    writeLines(paste(c(gene$chrom, "GeneTrackR_demo", "gene", gene$gene_start, gene$gene_end, ".", gene$strand, ".", gene_attr), collapse = "\t"), con)
    tx_set <- transcripts[transcripts$gene_id == gene$gene_id, , drop = FALSE]
    tx_set <- tx_set[order(tx_set$transcript_id), , drop = FALSE]
    for (j in seq_len(nrow(tx_set))) {
      tx <- tx_set[j, , drop = FALSE]
      tx_attr <- paste0("ID=", tx$transcript_id, ";Parent=", tx$gene_id, ";Name=", tx$transcript_id, ";transcript_biotype=", tx$gene_type)
      writeLines(paste(c(tx$chrom, "GeneTrackR_demo", "mRNA", tx$tx_start, tx$tx_end, ".", tx$strand, ".", tx_attr), collapse = "\t"), con)
      sub <- build_subfeatures(tx)
      type_count <- integer()
      names(type_count) <- character()
      for (k in seq_len(nrow(sub))) {
        type <- as.character(sub$type[k])
        if (!type %in% names(type_count)) type_count[type] <- 0L
        type_count[type] <- type_count[type] + 1L
        feature_id <- paste(tx$transcript_id, type, type_count[type], sep = ".")
        attr <- paste0(
          "ID=", feature_id,
          ";Parent=", tx$transcript_id,
          ";gene_id=", tx$gene_id,
          ";transcript_id=", tx$transcript_id,
          ";exon_number=", sub$exon_number[k]
        )
        writeLines(paste(c(tx$chrom, "GeneTrackR_demo", type, sub$start[k], sub$end[k], ".", tx$strand, sub$phase[k], attr), collapse = "\t"), con)
      }
    }
  }
}

write_gtf <- function(genes, transcripts, file) {
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  for (i in seq_len(nrow(genes))) {
    gene <- genes[i, , drop = FALSE]
    gene_attr <- paste0('gene_id "', gene$gene_id, '"; gene_name "', gene$gene_id, '"; gene_biotype "', gene$gene_type, '";')
    writeLines(paste(c(gene$chrom, "GeneTrackR_demo", "gene", gene$gene_start, gene$gene_end, ".", gene$strand, ".", gene_attr), collapse = "\t"), con)
    tx_set <- transcripts[transcripts$gene_id == gene$gene_id, , drop = FALSE]
    tx_set <- tx_set[order(tx_set$transcript_id), , drop = FALSE]
    for (j in seq_len(nrow(tx_set))) {
      tx <- tx_set[j, , drop = FALSE]
      tx_attr <- paste0(
        'gene_id "', tx$gene_id, '"; transcript_id "', tx$transcript_id,
        '"; gene_name "', tx$gene_id, '"; gene_biotype "', tx$gene_type,
        '"; transcript_biotype "', tx$gene_type, '";'
      )
      writeLines(paste(c(tx$chrom, "GeneTrackR_demo", "transcript", tx$tx_start, tx$tx_end, ".", tx$strand, ".", tx_attr), collapse = "\t"), con)
      sub <- build_subfeatures(tx)
      for (k in seq_len(nrow(sub))) {
        attr <- paste0(
          'gene_id "', tx$gene_id, '"; transcript_id "', tx$transcript_id,
          '"; gene_name "', tx$gene_id, '"; exon_number "', sub$exon_number[k],
          '"; gene_biotype "', tx$gene_type, '";'
        )
        writeLines(paste(c(tx$chrom, "GeneTrackR_demo", sub$type[k], sub$start[k], sub$end[k], ".", tx$strand, sub$phase[k], attr), collapse = "\t"), con)
      }
    }
  }
}

pattern_gt <- function(pattern, sample_id, hap_group, within_group) {
  sample_index <- as.integer(sub("^S", "", sample_id))
  if (pattern == "p13") state <- hap_group %in% c(3L, 4L)
  else if (pattern == "p23") state <- hap_group %in% c(2L, 3L)
  else if (pattern == "p4") state <- hap_group == 4L
  else if (pattern == "p2") state <- hap_group == 2L
  else if (pattern == "p14") state <- hap_group %in% c(1L, 4L)
  else if (pattern == "even") state <- sample_index %% 2L == 0L
  else if (pattern == "odd") state <- sample_index %% 2L == 1L
  else if (pattern == "mod3") state <- sample_index %% 3L == 0L
  else if (pattern == "mod4") state <- sample_index %% 4L %in% c(1L, 2L)
  else if (pattern == "blocks3") state <- ((sample_index - 1L) %/% 3L) %% 2L == 1L
  else if (grepl("^ldgrad[0-9]{2}$", pattern)) {
    flip_n <- as.integer(sub("^ldgrad", "", pattern))
    stopifnot(!is.na(flip_n), flip_n >= 0L, flip_n <= 36L)
    state <- hap_group %in% c(3L, 4L)
    flip_rank <- if (sample_index <= 18L) {
      2L * sample_index - 1L
    } else {
      2L * (sample_index - 18L)
    }
    if (flip_rank <= flip_n) state <- !state
  } else if (pattern == "hetero_sparse") {
    if (sample_index %in% c(5L, 14L, 23L, 32L)) return("0/1")
    state <- hap_group %in% c(2L, 4L)
  } else if (pattern == "mixed_missing") {
    if (sample_index %in% c(3L, 12L, 21L, 30L)) return("./.")
    if (sample_index %in% c(5L, 14L, 23L, 32L)) return("0/1")
    state <- hap_group %in% c(3L, 4L)
  } else {
    stop("Unknown genotype pattern: ", pattern, call. = FALSE)
  }
  if (isTRUE(state)) "1/1" else "0/0"
}

write_vcf <- function(variants, samples, chromosomes, file) {
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("##fileformat=VCFv4.2", con)
  writeLines("##source=GeneTrackR_demo_v0.5.32", con)
  for (i in seq_len(nrow(chromosomes))) {
    writeLines(paste0("##contig=<ID=", chromosomes$chrom[i], ",length=", chromosomes$size[i], ">"), con)
  }
  writeLines('##INFO=<ID=ROLE,Number=1,Type=String,Description="Designed role in GeneTrackR demo data">', con)
  writeLines('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', con)
  writeLines(paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples$sample_id), collapse = "\t"), con)
  for (i in seq_len(nrow(variants))) {
    variant <- variants[i, , drop = FALSE]
    gt <- vapply(seq_len(nrow(samples)), function(j) {
      pattern_gt(variant$pattern, samples$sample_id[j], samples$hap_group[j], samples$within_group[j])
    }, character(1L))
    info <- paste0("ROLE=", gsub(";", ",", variant$role, fixed = TRUE))
    fields <- c(variant$chrom, variant$pos, variant$variant_id, variant$ref, variant$alt, 60, "PASS", info, "GT", gt)
    writeLines(paste(fields, collapse = "\t"), con)
  }
}

write_phenotype <- function(samples, file) {
  residual <- c(-1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, -0.6, 0.6)
  flowering_residual <- c(-2, -1, 0, 1, 2, -1, 0, 1, 0)
  seed_base <- c(`1` = 20, `2` = 25, `3` = 31, `4` = 23)
  height_base <- c(`1` = 100, `2` = 104, `3` = 110, `4` = 106)
  protein_base <- c(`1` = 38, `2` = 44, `3` = 44, `4` = 38)
  x <- samples
  x$seed_weight <- seed_base[as.character(x$hap_group)] + residual[x$within_group]
  x$protein_content <- protein_base[as.character(x$hap_group)] + residual[x$within_group] * 0.35
  x$plant_height <- height_base[as.character(x$hap_group)] + residual[x$within_group] * 1.5
  x$flowering_time <- 45 + flowering_residual[x$within_group]
  x$flower_color <- ifelse(x$within_group %in% c(1L, 3L, 5L, 7L, 9L), "Purple", "White")
  order_index <- c(19, 1, 28, 10, 20, 2, 29, 11, 21, 3, 30, 12, 22, 4, 31, 13, 23, 5, 32, 14, 24, 6, 33, 15, 25, 7, 34, 16, 26, 8, 35, 17, 27, 9, 36, 18)
  x <- x[match(sprintf("S%02d", order_index), x$sample_id), , drop = FALSE]
  out <- x[, c("sample_id", "seed_weight", "protein_content", "plant_height", "flowering_time", "flower_color")]
  utils::write.table(out, file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

get_cds_segments <- function(tx_row) {
  if (is.na(tx_row$cds_start) || is.na(tx_row$cds_end) ||
      !nzchar(as.character(tx_row$cds_start)) || !nzchar(as.character(tx_row$cds_end))) {
    return(data.frame(start = integer(), end = integer()))
  }
  exons <- get_exons(tx_row)
  cds_start <- as.integer(tx_row$cds_start)
  cds_end <- as.integer(tx_row$cds_end)
  out <- lapply(seq_len(nrow(exons)), function(i) {
    start <- max(exons$start[i], cds_start)
    end <- min(exons$end[i], cds_end)
    if (start > end) return(NULL)
    data.frame(start = start, end = end)
  })
  out <- Filter(Negate(is.null), out)
  if (length(out) == 0L) return(data.frame(start = integer(), end = integer()))
  do.call(rbind, out)
}

build_rnaseq_signal <- function(transcripts, signal_design, strand = c("+", "-")) {
  strand <- match.arg(strand)
  exon_factors <- c(1.00, 0.90, 1.10, 0.95)
  local_wave <- c(0.92, 1.00, 1.08, 0.96, 1.04)
  rows <- list()
  k <- 0L

  for (i in seq_len(nrow(signal_design))) {
    weight <- as.numeric(signal_design$rnaseq_weight[i])
    if (is.na(weight) || weight <= 0) next
    tx <- transcripts[transcripts$transcript_id == signal_design$transcript_id[i], , drop = FALSE]
    if (tx$strand[1L] != strand) next
    if (nrow(tx) != 1L) {
      stop("Signal design transcript was not found uniquely: ", signal_design$transcript_id[i], call. = FALSE)
    }
    exons <- get_exons(tx)
    for (j in seq_len(nrow(exons))) {
      positions <- seq.int(exons$start[j], exons$end[j])
      block_index <- (positions - exons$start[j]) %/% 100L
      values <- weight * exon_factors[(j - 1L) %% length(exon_factors) + 1L] *
        local_wave[block_index %% length(local_wave) + 1L]
      k <- k + 1L
      rows[[k]] <- data.frame(
        chrom = tx$chrom,
        pos = positions,
        value = values,
        stringsAsFactors = FALSE
      )
    }
  }

  signal <- do.call(rbind, rows)
  signal <- stats::aggregate(value ~ chrom + pos, data = signal, FUN = sum)
  signal$value <- round(signal$value, 3)
  signal[order(signal$chrom, signal$pos), , drop = FALSE]
}

compress_rnaseq_signal <- function(signal) {
  if (nrow(signal) == 0L) {
    return(data.frame(chrom = character(), start = integer(), end = integer(), value = numeric()))
  }
  signal <- signal[order(signal$chrom, signal$pos), , drop = FALSE]
  new_run <- c(
    TRUE,
    signal$chrom[-1L] != signal$chrom[-nrow(signal)] |
      signal$pos[-1L] != signal$pos[-nrow(signal)] + 1L |
      signal$value[-1L] != signal$value[-nrow(signal)]
  )
  run_id <- cumsum(new_run)
  groups <- split(seq_len(nrow(signal)), run_id)
  out <- lapply(groups, function(idx) {
    data.frame(
      chrom = signal$chrom[idx[1L]],
      start = signal$pos[idx[1L]] - 1L,
      end = signal$pos[idx[length(idx)]],
      value = signal$value[idx[1L]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

write_rnaseq_signal <- function(file, transcripts, signal_design, strand = c("+", "-")) {
  strand <- match.arg(strand)
  signal <- build_rnaseq_signal(transcripts, signal_design, strand = strand)
  bedgraph <- compress_rnaseq_signal(signal)
  utils::write.table(
    bedgraph,
    file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

build_riboseq_signal <- function(transcripts, signal_design, strand = c("+", "-")) {
  strand <- match.arg(strand)

  # Generate moderately dense integer P-site-like counts with heterogeneous
  # heights within every frame. Total counts are calibrated to approximately
  # 80:10:10 for frame0:frame1:frame2. Frame 0 is broadly occupied and carries
  # the dominant peaks, whereas frame 1 and frame 2 occur at distinct subsets
  # of codons with lower, irregular integer counts.
  frame0_wave <- c(
    0.55, 0.80, 1.10, 1.40, 0.70, 1.55, 0.95,
    1.25, 0.60, 1.35, 1.00, 1.50, 0.85
  )
  frame0_jitter <- c(0.90, 1.12, 0.82, 1.25, 1.00, 1.18, 0.94, 1.08, 0.86)
  rows <- list()
  k <- 0L

  allocate_integer_counts <- function(weights, target_total) {
    weights <- as.numeric(weights)
    weights[!is.finite(weights) | weights < 0] <- 0
    target_total <- as.integer(round(target_total))
    out <- integer(length(weights))
    if (target_total <= 0L || sum(weights) <= 0) {
      return(out)
    }

    scaled <- weights / sum(weights) * target_total
    out <- as.integer(floor(scaled))
    remainder <- target_total - sum(out)
    if (remainder > 0L) {
      fractional <- scaled - out
      order_index <- order(
        -fractional,
        -weights,
        seq_along(weights)
      )
      selected <- order_index[seq_len(remainder)]
      out[selected] <- out[selected] + 1L
    }
    out
  }

  for (i in seq_len(nrow(signal_design))) {
    weight <- as.numeric(signal_design$riboseq_weight[i])
    if (is.na(weight) || weight <= 0) next

    tx <- transcripts[transcripts$transcript_id == signal_design$transcript_id[i], , drop = FALSE]
    if (tx$strand[1L] != strand) next
    if (nrow(tx) != 1L) {
      stop("Signal design transcript was not found uniquely: ", signal_design$transcript_id[i], call. = FALSE)
    }

    cds <- get_cds_segments(tx)
    if (nrow(cds) == 0L) next

    if (tx$strand == "+") {
      cds <- cds[order(cds$start), , drop = FALSE]
      positions <- unlist(
        lapply(seq_len(nrow(cds)), function(j) seq.int(cds$start[j], cds$end[j])),
        use.names = FALSE
      )
    } else {
      cds <- cds[order(cds$start, decreasing = TRUE), , drop = FALSE]
      positions <- unlist(
        lapply(seq_len(nrow(cds)), function(j) seq.int(cds$end[j], cds$start[j])),
        use.names = FALSE
      )
    }

    if (length(positions) %% 3L != 0L) {
      stop(
        "Primary coding transcript CDS length is not divisible by 3: ",
        signal_design$transcript_id[i],
        call. = FALSE
      )
    }

    n_codons <- length(positions) %/% 3L
    codon_index <- seq_len(n_codons) - 1L

    # Frame 0 is present at most codons but has strongly heterogeneous heights.
    # The occupancy rule is deterministic so regenerated demo data are stable.
    frame0_active <- ((codon_index * 19L + 7L) %% 31L) < 28L
    frame0_active[1L] <- TRUE
    frame0_active[n_codons] <- TRUE
    frame0_raw <- weight *
      frame0_wave[codon_index %% length(frame0_wave) + 1L] *
      frame0_jitter[codon_index %% length(frame0_jitter) + 1L]
    frame0 <- integer(n_codons)
    frame0[frame0_active] <- pmax(
      1L,
      as.integer(round(frame0_raw[frame0_active]))
    )

    internal_frame0 <- frame0
    internal_frame0[c(1L, n_codons)] <- 0L
    internal_mean <- mean(internal_frame0[internal_frame0 > 0L])
    if (!is.finite(internal_mean) || internal_mean <= 0) {
      internal_mean <- weight
    }
    boundary_count <- max(1L, as.integer(round(2 * internal_mean)))
    frame0[1L] <- boundary_count
    frame0[n_codons] <- boundary_count

    # Each off-frame receives one eighth of the frame-0 total, giving an
    # approximately 80:10:10 total-count ratio. Frame 1 and frame 2 use
    # different active codons and different pseudo-irregular weight patterns,
    # so their bars are neither synchronized nor nearly equal in height.
    off_frame_total <- as.integer(round(sum(frame0) / 8))

    frame1_active <- ((codon_index * 13L + 5L) %% 29L) < 13L
    frame2_active <- ((codon_index * 17L + 9L) %% 31L) < 14L
    frame1_active[c(1L, n_codons)] <- FALSE
    frame2_active[c(1L, n_codons)] <- FALSE

    frame1_weights <- frame1_active *
      (1 + ((codon_index * 7L + 3L) %% 9L)) *
      (0.75 + 0.08 * ((codon_index * 5L + 2L) %% 7L))
    frame2_weights <- frame2_active *
      (1 + ((codon_index * 11L + 4L) %% 8L)) *
      (0.80 + 0.07 * ((codon_index * 3L + 1L) %% 8L))

    frame1 <- allocate_integer_counts(frame1_weights, off_frame_total)
    frame2 <- allocate_integer_counts(frame2_weights, off_frame_total)

    values <- numeric(length(positions))
    values[seq.int(1L, length(values), by = 3L)] <- frame0
    values[seq.int(2L, length(values), by = 3L)] <- frame1
    values[seq.int(3L, length(values), by = 3L)] <- frame2

    keep <- values > 0
    if (!any(keep)) next

    k <- k + 1L
    rows[[k]] <- data.frame(
      chrom = tx$chrom,
      pos = as.integer(positions[keep]),
      value = as.numeric(values[keep]),
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0L) {
    return(data.frame(chrom = character(), pos = integer(), value = numeric()))
  }

  signal <- do.call(rbind, rows)
  signal <- stats::aggregate(value ~ chrom + pos, data = signal, FUN = sum)
  signal$value <- as.numeric(round(signal$value))
  signal <- signal[signal$value > 0, , drop = FALSE]
  signal[order(signal$chrom, signal$pos), , drop = FALSE]
}

write_riboseq_signal <- function(file, transcripts, signal_design, strand = c("+", "-")) {
  strand <- match.arg(strand)
  signal <- build_riboseq_signal(transcripts, signal_design, strand = strand)
  bedgraph <- data.frame(
    chrom = signal$chrom,
    start = signal$pos - 1L,
    end = signal$pos,
    value = signal$value,
    stringsAsFactors = FALSE
  )
  utils::write.table(
    bedgraph,
    file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

script_dir <- get_script_dir()
model_dir <- file.path(script_dir, "demo_model")
out_dir <- normalizePath(file.path(script_dir, "..", "extdata"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

old_signal_files <- c(
  "gtr_demo_ctrl_plus.bedgraph",
  "gtr_demo_ctrl_minus.bedgraph",
  "gtr_demo_treat_plus.bedgraph",
  "gtr_demo_treat_minus.bedgraph",
  "gtr_demo_rnaseq.bedgraph",
  "gtr_demo_riboseq.bedgraph",
  "gtr_demo_rnaseq_plus.bedgraph",
  "gtr_demo_rnaseq_minus.bedgraph",
  "gtr_demo_riboseq_plus.bedgraph",
  "gtr_demo_riboseq_minus.bedgraph"
)
unlink(file.path(out_dir, old_signal_files), force = TRUE)

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
unlink(file.path(out_dir, legacy_example_files), force = TRUE)

chromosomes <- read_tsv(file.path(model_dir, "chromosomes.tsv"))
transcripts <- read_tsv(file.path(model_dir, "transcripts.tsv"))
variants <- read_tsv(file.path(model_dir, "variants.tsv"))
samples <- read_tsv(file.path(model_dir, "samples.tsv"))
features <- read_tsv(file.path(model_dir, "features.tsv"))
signal_design <- read_tsv(file.path(model_dir, "signal_design.tsv"))

transcripts$tx_start <- as.integer(transcripts$tx_start)
transcripts$tx_end <- as.integer(transcripts$tx_end)
transcripts$cds_start <- suppressWarnings(as.integer(transcripts$cds_start))
transcripts$cds_end <- suppressWarnings(as.integer(transcripts$cds_end))
variants$pos <- as.integer(variants$pos)
samples$hap_group <- as.integer(samples$hap_group)
samples$within_group <- as.integer(samples$within_group)
features$start <- as.integer(features$start)
features$end <- as.integer(features$end)
chromosomes$size <- as.integer(chromosomes$size)
signal_design$rnaseq_weight <- as.numeric(signal_design$rnaseq_weight)
signal_design$riboseq_weight <- as.numeric(signal_design$riboseq_weight)
genes <- build_gene_table(transcripts)

write_genepred_ext(transcripts, file.path(out_dir, "gtr_demo.genePredExt"))
write_gff3(genes, transcripts, file.path(out_dir, "gtr_demo.gff3"))
write_gtf(genes, transcripts, file.path(out_dir, "gtr_demo.gtf"))

bed <- data.frame(
  chrom = features$chrom,
  start = features$start - 1L,
  end = features$end,
  name = paste(features$feature_id, features$type, sep = "|"),
  score = 0,
  strand = features$strand,
  stringsAsFactors = FALSE
)
utils::write.table(bed, file.path(out_dir, "gtr_demo_features.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
utils::write.table(chromosomes, file.path(out_dir, "gtr_demo.chrom.sizes"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write_vcf(variants, samples, chromosomes, file.path(out_dir, "gtr_demo_variants.vcf"))
write_phenotype(samples, file.path(out_dir, "gtr_demo_pheno.tsv"))
write_rnaseq_signal(file.path(out_dir, "gtr_demo_rnaseq_plus.bedgraph"), transcripts, signal_design, strand = "+")
write_rnaseq_signal(file.path(out_dir, "gtr_demo_rnaseq_minus.bedgraph"), transcripts, signal_design, strand = "-")
write_riboseq_signal(file.path(out_dir, "gtr_demo_riboseq_plus.bedgraph"), transcripts, signal_design, strand = "+")
write_riboseq_signal(file.path(out_dir, "gtr_demo_riboseq_minus.bedgraph"), transcripts, signal_design, strand = "-")

message("[GeneTrackR demo] Generated files in: ", out_dir)
message("[GeneTrackR demo] Genes: ", nrow(genes), "; transcripts: ", nrow(transcripts), "; samples: ", nrow(samples), "; variants: ", nrow(variants), ".")
