#!/usr/bin/env Rscript

# Author: Rensc
# Date: 2026-11-21
# Version: dev001
# Function: Compare TxDb relationships and indexed VCF readers on fixed fixtures.
# Input: Installed GeneTrackR, optional readers, a local report directory.
# Output: Key-level mapping, count, timing and peak-memory CSV reports.

membership_rows <- function(kind, transcript, chrom, start, end, strand) {
  data.frame(kind = rep(kind, length(transcript)), transcript = as.character(transcript),
    chrom = as.character(chrom), start = as.integer(start), end = as.integer(end),
    strand = as.character(strand))
}

range_membership <- function(kind, transcript, ranges) {
  membership_rows(kind, transcript, GenomicRanges::seqnames(ranges),
    IRanges::start(ranges), IRanges::end(ranges), GenomicRanges::strand(ranges))
}

assess_txdb <- function(file, output) {
  feature <- GeneTrackR::read_gff(file, verbose = FALSE, progress = FALSE)$data
  txdb <- txdbmaker::makeTxDbFromGFF(file, format = "gff3")
  kinds <- c(transcript = "mRNA", exon = "exon", cds = "CDS")
  gff_membership <- function(use_all_parents) {
    do.call(rbind, lapply(names(kinds), function(kind) {
      rows <- feature[feature$type == kinds[[kind]], ]
      parents <- if (kind == "transcript" || !use_all_parents) {
        as.list(rows$transcript_id)
      } else {
        rows$parent_ids
      }
      index <- rep(seq_len(nrow(rows)), lengths(parents))
      membership_rows(kind, unlist(parents, use.names = FALSE), rows$chrom[index],
        rows$start[index], rows$end[index], rows$strand[index])
    }))
  }
  grouped_membership <- function(kind, groups) {
    range_membership(kind, rep(names(groups), lengths(groups)), unlist(groups, use.names = FALSE))
  }
  transcripts <- GenomicFeatures::transcripts(txdb, columns = "tx_name")
  standard <- rbind(
    range_membership("transcript", transcripts$tx_name, transcripts),
    grouped_membership("exon", GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)),
    grouped_membership("cds", GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE))
  )
  declared <- gff_membership(TRUE)
  primary <- gff_membership(FALSE)
  key <- function(x) do.call(paste, c(x, sep = "\t"))
  raw_standard <- standard
  # Compare TxDb display names to GFF3 IDs through the fixture's Name values.
  tx_rows <- feature[feature$type == kinds[["transcript"]], ]
  name_map <- data.frame(gff3_id = tx_rows$feature_id, txdb_name = tx_rows$name)
  stopifnot(!anyNA(name_map), !anyDuplicated(name_map$gff3_id), !anyDuplicated(name_map$txdb_name))
  mapped <- match(standard$transcript, name_map$txdb_name)
  stopifnot(!anyNA(mapped))
  standard$transcript <- name_map$gff3_id[mapped]
  for (rows in list(declared, primary, standard)) {
    stopifnot(!anyNA(rows), !anyDuplicated(key(rows)))
  }
  all_rows <- unique(rbind(declared, primary, standard))
  all_rows$in_gff3 <- key(all_rows) %in% key(declared)
  all_rows$in_primary <- key(all_rows[names(declared)]) %in% key(primary)
  all_rows$in_txdb <- key(all_rows[names(declared)]) %in% key(standard)
  counts <- do.call(rbind, lapply(names(kinds), function(kind) {
    compared <- all_rows[all_rows$kind == kind, ]
    data.frame(kind = kind, source_rows = sum(feature$type == kinds[[kind]]),
      declared_links = sum(compared$in_gff3), primary_links = sum(compared$in_primary),
      txdb_links = sum(compared$in_txdb),
      matches_before_name_mapping = sum(key(compared[compared$in_gff3, names(declared)]) %in% key(raw_standard)),
      matched_links = sum(compared$in_gff3 & compared$in_txdb),
      gff3_only = sum(compared$in_gff3 & !compared$in_txdb),
      txdb_only = sum(compared$in_txdb & !compared$in_gff3),
      declared_missing_in_primary = sum(compared$in_gff3 & !compared$in_primary))
  }))
  utils::write.csv(all_rows, file.path(output, "txdb-rows.csv"), row.names = FALSE)
  utils::write.csv(name_map, file.path(output, "txdb-name-map.csv"), row.names = FALSE)
  utils::write.csv(counts, file.path(output, "txdb-counts.csv"), row.names = FALSE)
  utils::write.csv(as.data.frame(table(feature$type)), file.path(output, "gff3-records.csv"),
    row.names = FALSE)
  print(counts, row.names = FALSE)
}

peak_working_set_mib <- function() {
  if (.Platform$OS.type == "windows") {
    query <- sprintf("(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid())
    bytes <- system2("powershell.exe", c("-NoProfile", "-NonInteractive", "-Command", shQuote(query)),
      stdout = TRUE)
    stopifnot(length(bytes) == 1L, !is.na(as.numeric(bytes)))
    return(as.numeric(bytes) / 1024^2)
  }
  if (file.exists("/proc/self/status")) {
    status <- readLines("/proc/self/status")
    high_water <- grep("^VmHWM:", status, value = TRUE)
    if (length(high_water)) return(as.numeric(gsub("[^0-9]", "", high_water)) / 1024)
  }
  NA_real_
}

indexed_vcf_fixture <- function(output, records = 50000L) {
  samples <- sprintf("S%02d", c(2L, 1L, seq.int(3L, 8L)))
  positions <- as.integer(seq.int(1L, by = 10L, length.out = records))
  reference <- ifelse(positions == 99991L, strrep("A", 15L), "A")
  alternate <- ifelse(positions == 99991L, "A", "C")
  sites <- data.frame(chrom = "chr1", pos = positions, id = paste0("v", positions),
    ref = reference, alt = alternate, depth = positions %% 19L + 1L)
  gt <- outer(seq_along(positions), seq_along(samples), function(row, col) {
    c("0|0", "0|1", "1|1")[(row + col) %% 3L + 1L]
  })
  region <- c(start = 100001L, end = 100500L)
  selected <- positions >= region[["start"]] & positions <= region[["end"]]
  input <- file.path(output, "regional-query.vcf")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=1000000>",
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples),
      collapse = "\t")), input)
  rows <- cbind(sites[setdiff(names(sites), "depth")], 60, "PASS", paste0("DP=", sites$depth), "GT", gt)
  utils::write.table(rows, input, append = TRUE, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE)
  indexed <- Rsamtools::bgzip(input, dest = paste0(input, ".gz"))
  Rsamtools::indexTabix(indexed, format = "vcf")
  list(input = indexed, records = records, samples = samples, region = region,
    sites = sites[selected, ], gt = gt[selected, , drop = FALSE])
}

assess_indexed_worker <- function(reader, output) {
  fixture <- readRDS(file.path(output, "indexed-fixture.rds"))
  initial <- loadedNamespaces()
  package <- if (reader == "native") "GeneTrackR" else "VariantAnnotation"
  startup <- system.time({
    loadNamespace(package)
    loadNamespace("Rsamtools")
  })[["elapsed"]]
  loaded <- length(setdiff(loadedNamespaces(), initial))
  start <- fixture$region[["start"]]
  end <- fixture$region[["end"]]
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start, end))
  query <- if (reader == "native") {
    function() GeneTrackR::read_vcf(fixture$input, chrom = "chr1", start = start, end = end,
      verbose = FALSE, progress = FALSE)
  } else {
    param <- VariantAnnotation::ScanVcfParam(which = region)
    function() VariantAnnotation::readVcf(fixture$input, genome = NA_character_, param = param)
  }
  timings <- numeric(3L)
  for (i in seq_along(timings)) {
    gc()
    timings[[i]] <- system.time(result <- query())[["elapsed"]]
  }
  peak <- peak_working_set_mib()
  if (reader == "native") {
    table <- result$data
    samples <- result$meta$sample_names
    sites <- data.frame(chrom = table$chrom, pos = table$pos, id = table$variant_id,
      ref = table$ref, alt = table$alt, depth = as.integer(sub("^DP=", "", table$info)))
    gt <- as.matrix(table[, samples, with = FALSE])
  } else {
    ranges <- GenomicRanges::granges(result)
    fixed <- VariantAnnotation::fixed(result)
    samples <- colnames(result)
    sites <- data.frame(chrom = as.character(GenomicRanges::seqnames(ranges)),
      pos = IRanges::start(ranges), id = rownames(result), ref = as.character(fixed$REF),
      alt = vapply(as.list(fixed$ALT), paste, character(1L), collapse = ","),
      depth = VariantAnnotation::info(result)$DP)
    gt <- VariantAnnotation::geno(result)$GT
  }
  keep <- sites$chrom == "chr1" & sites$pos >= start & sites$pos <= end
  selected <- sites[keep, ]
  expected <- fixture$sites
  stopifnot(!anyDuplicated(selected$id), !anyDuplicated(expected$id))
  index <- match(expected$id, selected$id)
  matched <- !is.na(index)
  field_matches <- vapply(names(expected), function(field) {
    identical(expected[[field]][matched], selected[[field]][index[matched]])
  }, logical(1L))
  same_samples <- identical(samples, fixture$samples)
  same_gt <- identical(unname(gt[keep, , drop = FALSE][index[matched], , drop = FALSE]),
    unname(fixture$gt[matched, , drop = FALSE]))
  counts <- data.frame(reader = reader, version = as.character(getNamespaceVersion(package)),
    input_sites = fixture$records, input_samples = length(fixture$samples),
    raw_returned = nrow(sites), outside_pos = sum(!keep), expected_pos = nrow(expected),
    matched_pos = sum(matched), expected_only = sum(!matched),
    reader_only = sum(!selected$id %in% expected$id), sites_match = all(field_matches),
    samples_match = same_samples, genotypes_match = same_gt, namespaces_loaded = loaded,
    startup_seconds = startup, median_read_seconds = stats::median(timings), peak_rss_mib = peak)
  utils::write.csv(counts, file.path(output, paste0("indexed-", reader, ".csv")), row.names = FALSE)
  utils::write.csv(data.frame(reader = reader, elapsed_seconds = timings),
    file.path(output, paste0("timings-", reader, ".csv")), row.names = FALSE)
  stopifnot(counts$expected_only == 0L, counts$reader_only == 0L,
    all(field_matches), same_samples, same_gt)
}

assess_standard_backends <- function(output) {
  if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE))) {
    stop("Choose an empty report directory.", call. = FALSE)
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, winslash = "/", mustWork = TRUE)
  gff <- "tools/io-evaluation/fixtures/annotation.gff3"
  assess_txdb(gff, output)
  fixture <- indexed_vcf_fixture(output)
  saveRDS(fixture, file.path(output, "indexed-fixture.rds"))
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  rscript <- file.path(R.home("bin"), "Rscript")
  for (reader in c("native", "standard")) {
    status <- system2(rscript, c("--vanilla", shQuote(script), "--worker", reader, shQuote(output)))
    if (status != 0L) stop(sprintf("The %s indexed reader assessment failed.", reader), call. = FALSE)
  }
  packages <- c("GeneTrackR", "Rsamtools", "txdbmaker", "GenomicFeatures", "VariantAnnotation")
  versions <- vapply(packages, function(package) {
    loadNamespace(package)
    as.character(getNamespaceVersion(package))
  }, character(1L))
  details <- c(Platform = unname(Sys.info()[["sysname"]]), R = as.character(getRversion()), versions,
    GFF3_MD5 = unname(tools::md5sum(gff)), VCF_MD5 = unname(tools::md5sum(sub("[.]gz$", "", fixture$input))))
  write.dcf(as.data.frame(as.list(details)), file.path(output, "metadata.dcf"))
  evidence <- do.call(rbind, lapply(c("native", "standard"), function(reader) {
    utils::read.csv(file.path(output, paste0("indexed-", reader, ".csv")))
  }))
  utils::write.csv(evidence, file.path(output, "indexed-vcf.csv"), row.names = FALSE)
  print(evidence, row.names = FALSE)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 3L && identical(args[[1L]], "--worker") && args[[2L]] %in% c("native", "standard")) {
    assess_indexed_worker(args[[2L]], args[[3L]])
  } else if (length(args) == 1L) {
    assess_standard_backends(args[[1L]])
  } else {
    stop("Usage: Rscript tools/assess_standard_backends.R OUTPUT_DIR", call. = FALSE)
  }
}
