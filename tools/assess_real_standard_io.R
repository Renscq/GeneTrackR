# Author: Rensc
# Date: 2026-11-21
# Version: dev001
# Function: Compare native and standard readers on pinned public data.
# Input: A directory with the chr22 VCF, its index, and the GENCODE v19 GFF3.
# Output: Reader timings, peak RSS, and row, sample, and field comparisons.

real_io_inputs <- c(
  "chr22.vcf.gz" = "46e78a67f51161387c4c13aeefe80367",
  "chr22.vcf.gz.tbi" = "540ee7e114a2c3c9f9730610940ec47f",
  "gencode.v19.annotation.gff3.gz" = "89ca466f30dd01539ac48d1e8b7964ba"
)
real_io_region <- c(chrom = "22", start = "50443100", end = "51000000")
real_io_readers <- c(
  vcf_native = "GeneTrackR", vcf_standard = "VariantAnnotation",
  gff3_native = "GeneTrackR", gff3_standard = "rtracklayer"
)

real_io_granges <- function() {
  GenomicRanges::GRanges(
    real_io_region[["chrom"]],
    IRanges::IRanges(
      as.integer(real_io_region[["start"]]), as.integer(real_io_region[["end"]])
    )
  )
}

extract_chr22_gff3 <- function(path, output) {
  input <- gzfile(path, "rt")
  on.exit(close(input), add = TRUE)
  target <- file(output, "wb")
  on.exit(close(target), add = TRUE)
  writeLines("##gff-version 3", target)
  counts <- c(source_rows = 0L, selected_rows = 0L)
  repeat {
    lines <- readLines(input, n = 10000L)
    if (!length(lines)) break
    keep <- startsWith(lines, "chr22\t")
    counts <- counts + c(sum(nzchar(lines) & !startsWith(lines, "#")), sum(keep))
    writeLines(lines[keep], target, useBytes = TRUE)
  }
  counts
}

real_io_worker <- function(reader, output) {
  inputs <- readRDS(file.path(output, "inputs.rds"))
  package <- real_io_readers[[reader]]
  startup <- system.time({
    loadNamespace(package)
    if (startsWith(reader, "vcf_")) loadNamespace("Rsamtools")
  })[["elapsed"]]
  query <- switch(reader,
    vcf_native = function() GeneTrackR::read_vcf(
      inputs$vcf, chrom = real_io_region[["chrom"]],
      start = as.integer(real_io_region[["start"]]),
      end = as.integer(real_io_region[["end"]]), verbose = FALSE, progress = FALSE
    ),
    vcf_standard = function() VariantAnnotation::readVcf(
      Rsamtools::TabixFile(inputs$vcf), "hg19",
      param = VariantAnnotation::ScanVcfParam(which = real_io_granges())
    ),
    gff3_native = function() GeneTrackR::read_gff(
      inputs$gff, verbose = FALSE, progress = FALSE
    ),
    gff3_standard = function() io_with_local_path(inputs$gff, function(path) {
      rtracklayer::import(path, format = "gff3")
    })
  )
  measured <- measure_io_reader(query)
  if (measured$status != "ok") stop(measured$error, call. = FALSE)
  saveRDS(measured$value, file.path(output, paste0(reader, ".rds")), compress = FALSE)
  data.frame(
    reader = reader, package = package,
    version = as.character(getNamespaceVersion(package)),
    startup_seconds = startup, measured$timing, peak_rss_mib = peak_working_set_mib(),
    warnings = measured$warnings
  ) |>
    utils::write.csv(file.path(output, paste0(reader, "-timing.csv")), row.names = FALSE)
}

run_real_io_reader <- function(reader, output, script) {
  rscript <- file.path(R.home("bin"), "Rscript")
  args <- c("--vanilla", shQuote(script), "--worker", reader, shQuote(output))
  log <- file.path(output, paste0(reader, ".log"))
  err <- file.path(output, paste0(reader, ".stderr"))
  mac <- identical(unname(Sys.info()[["sysname"]]), "Darwin")
  status <- if (mac) {
    system2("/usr/bin/time", c("-l", shQuote(rscript), args),
      stdout = log, stderr = err, env = "LC_ALL=C")
  } else {
    system2(rscript, args, stdout = log, stderr = err)
  }
  if (status != 0L) {
    stop(sprintf("Reader %s failed. See %s.", reader, err), call. = FALSE)
  }
  timings <- utils::read.csv(file.path(output, paste0(reader, "-timing.csv")))
  if (mac) {
    peak <- grep("^[[:space:]]*[0-9]+[[:space:]]+maximum resident set size",
      readLines(err, warn = FALSE), value = TRUE)
    timings$peak_rss_mib <- if (length(peak) == 1L) {
      as.numeric(sub("^\\s*([0-9]+).*", "\\1", peak)) / 1024^2
    } else NA_real_
  }
  timings
}

# Occurrences keep repeated identities distinct, in each reader's input order.
real_io_key <- function(table, fields) {
  columns <- lapply(table[fields], function(x) ifelse(is.na(x), "0", paste0("1", x)))
  key <- do.call(paste, c(unname(columns), sep = "\t"))
  paste(key, data.table::rowid(key), sep = "\t")
}

compare_real_rows <- function(source, observed, keys, fields, label, output) {
  source$key <- real_io_key(source, keys)
  observed$key <- real_io_key(observed, keys)
  result <- compare_io_tables(source, observed, "key", fields, tolerance = 1e-7)
  for (name in c("differences", "native_only", "standard_only")) {
    if (nrow(result[[name]])) {
      utils::write.csv(result[[name]], file.path(output, paste0(label, "-", name, ".csv")),
        row.names = FALSE)
    }
  }
  result <- result$summary
  names(result)[names(result) %in% c("native_rows", "standard_rows",
    "native_only_rows", "standard_only_rows")] <- c(
    "source_rows", "reader_rows", "source_only_rows", "reader_only_rows"
  )
  data.frame(comparison = label, result, source_order = identical(source$key, observed$key))
}

real_io_vcf_text <- function(lines, columns) {
  data.table::fread(text = paste(lines, collapse = "\n"), sep = "\t",
    header = FALSE, col.names = columns, colClasses = "character",
    na.strings = NULL, quote = "", data.table = FALSE, showProgress = FALSE)
}

vcf_text_values <- function(table, samples) {
  out <- table[c("chrom", "pos", "ref", "alt", "id", "qual", "filter")]
  out$pos <- as.integer(out$pos)
  out$qual <- as.numeric(io_missing(out$qual))
  out$filter <- io_missing(out$filter)
  out$id <- io_missing(out$id)
  out$AF <- as.numeric(io_info_values(table$info, "AF"))
  out$END <- as.integer(io_info_values(table$info, "END"))
  if (!all(table$FORMAT == "GT:DS:GL")) {
    stop("Unexpected FORMAT order in the pinned VCF.", call. = FALSE)
  }
  token <- function(parts, index) {
    vapply(parts, function(values) values[index], character(1L))
  }
  for (sample in samples) {
    parts <- strsplit(table[[sample]], ":", fixed = TRUE)
    out[[paste("GT", sample)]] <- io_missing(token(parts, 1L))
    out[[paste("DS", sample)]] <- as.numeric(io_missing(token(parts, 2L)))
    gl <- strsplit(token(parts, 3L), ",", fixed = TRUE)
    if (any(lengths(gl) > 3L)) stop("Unexpected VCF GL cardinality.", call. = FALSE)
    for (i in seq_len(3L)) {
      out[[paste0("GL", i, " ", sample)]] <- as.numeric(io_missing(token(gl, i)))
    }
  }
  out
}

vcf_standard_values <- function(object, samples) {
  ranges <- SummarizedExperiment::rowRanges(object)
  info <- VariantAnnotation::info(object)
  out <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(ranges)),
    pos = IRanges::start(ranges),
    ref = as.character(VariantAnnotation::ref(object)),
    alt = vapply(as.list(VariantAnnotation::alt(object)), paste, character(1L), collapse = ","),
    id = io_missing(names(ranges)), qual = VariantAnnotation::qual(object),
    filter = io_missing(VariantAnnotation::filt(object)), AF = info$AF, END = info$END
  )
  gt <- VariantAnnotation::geno(object)
  for (sample in samples) {
    out[[paste("GT", sample)]] <- io_missing(gt$GT[, sample])
    out[[paste("DS", sample)]] <- gt$DS[, sample]
    gl <- gt$GL[, sample]
    if (any(lengths(gl) > 3L)) stop("Unexpected standard GL cardinality.", call. = FALSE)
    for (i in seq_len(3L)) {
      out[[paste0("GL", i, " ", sample)]] <- vapply(gl, function(x) x[i], numeric(1L))
    }
  }
  out
}

compare_real_vcf <- function(path, native, standard, output) {
  raw_lines <- readLines(gzfile(path), warn = FALSE)
  columns <- strsplit(grep("^#CHROM\t", raw_lines, value = TRUE), "\t", fixed = TRUE)[[1L]]
  samples <- columns[seq.int(10L, length(columns))]
  columns[seq_len(8L)] <- c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info")
  raw_lines <- raw_lines[!startsWith(raw_lines, "#") & nzchar(raw_lines)]
  raw_vcf <- real_io_vcf_text(raw_lines, columns)
  tabix <- Rsamtools::scanTabix(Rsamtools::TabixFile(path), param = real_io_granges())[[1L]]
  chrom <- real_io_region[["chrom"]]
  start <- as.integer(real_io_region[["start"]])
  end <- as.integer(real_io_region[["end"]])
  inside <- function(table) {
    table$chrom == chrom & table$pos >= start & table$pos <= end
  }
  keys <- c("chrom", "pos", "ref", "alt")
  raw_pos <- as.integer(raw_vcf$pos)
  raw_end <- pmax(raw_pos + nchar(raw_vcf$ref) - 1L,
    as.integer(io_info_values(raw_vcf$info, "END")), na.rm = TRUE)
  overlap <- raw_vcf$chrom == chrom & raw_pos <= end & raw_end >= start
  source <- vcf_text_values(raw_vcf[overlap, ], samples)
  source_pos <- source[inside(source), , drop = FALSE]
  indexed <- real_io_vcf_text(tabix, columns)
  native_names <- native$meta$sample_names
  standard_names <- colnames(standard)
  same_samples <- c(native = identical(native_names, samples),
    standard = identical(standard_names, samples))
  utils::write.csv(data.frame(source = samples,
    native = match(samples, native_names), standard = match(samples, standard_names)),
    file.path(output, "samples.csv"), row.names = FALSE)
  if (anyDuplicated(samples) || !setequal(samples, native_names) ||
      !setequal(samples, standard_names)) {
    stop("VCF sample names differ from the source header.", call. = FALSE)
  }

  native_table <- as.data.frame(native$data)
  names(native_table)[names(native_table) == "variant_id"] <- "id"
  observed <- list(
    vcf_tabix = vcf_text_values(indexed, samples),
    vcf_native = vcf_text_values(native_table, samples),
    vcf_standard = vcf_standard_values(standard, samples)
  )
  counts <- c(vcf_source_rows = nrow(raw_vcf), vcf_sample_count = length(samples),
    vcf_source_overlap_rows = nrow(source), vcf_pos_rows = nrow(source_pos),
    vcf_missing_source_id = sum(is.na(source_pos$id)),
    vcf_duplicate_present_ids = sum(duplicated(stats::na.omit(source_pos$id))),
    vcf_duplicate_positions = sum(duplicated(source_pos[c("chrom", "pos")])),
    vcf_duplicate_allele_keys = sum(duplicated(source_pos[keys])),
    vcf_pos_rows_with_multibase_ref = sum(nchar(source_pos$ref) > 1L),
    vcf_overlap_rows_with_info_end = sum(!is.na(source$END)))
  fields <- setdiff(names(source), keys)
  comparisons <- vector("list", length(observed))
  source_key <- real_io_key(source, keys)
  for (i in seq_along(observed)) {
    reader <- names(observed)[i]
    table <- observed[[reader]]
    counts <- c(counts, stats::setNames(c(nrow(table), sum(!inside(table))),
      paste0(reader, c("_overlap_rows", "_outside_pos_rows"))))
    matched <- match(real_io_key(table, keys), source_key)
    missing_source_id <- !is.na(matched) & is.na(source$id[matched])
    counts <- c(counts, stats::setNames(
      sum(inside(table) & missing_source_id & !is.na(table$id)),
      paste0(reader, "_generated_ids")))
    table$id[missing_source_id] <- NA_character_
    comparisons[[i]] <- rbind(
      compare_real_rows(source, table, keys, fields, paste0(reader, "_overlap"), output),
      compare_real_rows(source_pos, table[inside(table), ], keys, fields,
        paste0(reader, "_pos"), output)
    )
  }
  list(comparisons = do.call(rbind, comparisons), counts = counts, same_samples = same_samples)
}

compare_real_gff <- function(path, native, standard, output) {
  source <- utils::read.delim(path, sep = "\t", comment.char = "#", quote = "",
    header = FALSE, col.names = c("chrom", "source", "type", "start", "end",
      "score", "strand", "phase", "attributes"),
    colClasses = c("character", "character", "character", "integer", "integer",
      "character", "character", "character", "character"))
  source$score <- as.numeric(io_missing(source$score))
  source$phase <- io_missing(source$phase)
  source$strand[source$strand == "."] <- "*"
  source$id <- io_info_values(source$attributes, "ID")
  source$parent <- io_info_values(source$attributes, "Parent")
  keys <- c("chrom", "start", "end", "strand", "type", "id", "parent")
  fields <- c("source", "score", "phase")
  observed <- list(
    gff3_native = io_annotation_table(native, "gff3", "native"),
    gff3_standard = io_annotation_table(standard, "gff3", "rtracklayer")
  )
  comparisons <- lapply(names(observed), function(reader) {
    compare_real_rows(source, observed[[reader]], keys, fields, reader, output)
  })
  utils::write.csv(as.data.frame(table(feature_type = source$type)),
    file.path(output, "gff3-feature-types.csv"), row.names = FALSE)
  relationship_counts <- vapply(c(list(gff3_source = source), observed), function(receiver) {
    sum(lengths(strsplit(stats::na.omit(receiver$parent), ",", fixed = TRUE)))
  }, integer(1L))
  list(comparisons = do.call(rbind, comparisons), counts = c(
    gff3_missing_ids = sum(is.na(source$id)),
    gff3_duplicate_present_ids = sum(duplicated(stats::na.omit(source$id))),
    stats::setNames(relationship_counts, paste0(names(relationship_counts), "_parent_links"))
  ))
}

assess_real_standard_io <- function(data_dir, output, script) {
  paths <- normalizePath(file.path(data_dir, names(real_io_inputs)), winslash = "/", mustWork = TRUE)
  checksums <- unname(tools::md5sum(paths))
  if (!identical(checksums, unname(real_io_inputs))) {
    stop("A benchmark input has an unexpected MD5 checksum.", call. = FALSE)
  }
  if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE))) {
    stop("Choose an empty report directory.", call. = FALSE)
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, winslash = "/", mustWork = TRUE)
  gff <- file.path(output, "gencode.v19.chr22.gff3")
  gff_counts <- extract_chr22_gff3(paths[3L], gff)
  saveRDS(list(vcf = paths[1L], gff = gff), file.path(output, "inputs.rds"))
  readers <- names(real_io_readers)
  timing <- do.call(rbind, lapply(readers, run_real_io_reader, output = output, script = script))
  utils::write.csv(timing, file.path(output, "timing.csv"), row.names = FALSE)
  objects <- stats::setNames(lapply(readers, function(reader) {
    readRDS(file.path(output, paste0(reader, ".rds")))
  }), readers)
  vcf <- compare_real_vcf(paths[1L], objects$vcf_native, objects$vcf_standard, output)
  annotation <- compare_real_gff(gff, objects$gff3_native, objects$gff3_standard, output)
  comparisons <- rbind(vcf$comparisons, annotation$comparisons)
  utils::write.csv(comparisons, file.path(output, "comparisons.csv"), row.names = FALSE)
  counts <- c(vcf$counts, stats::setNames(gff_counts, paste0("gff3_", names(gff_counts))),
    annotation$counts)
  utils::write.csv(data.frame(metric = names(counts), count = unname(counts)),
    file.path(output, "counts.csv"), row.names = FALSE)
  metadata <- c(
    Platform = Sys.info()[["sysname"]], Architecture = R.version$arch, R = R.version.string,
    Bioconductor = as.character(BiocManager::version()),
    Rsamtools = as.character(getNamespaceVersion("Rsamtools")),
    Region = sprintf("%s:%s-%s", real_io_region[1L], real_io_region[2L], real_io_region[3L]),
    VCF_MD5 = checksums[1L], TBI_MD5 = checksums[2L], GFF3_MD5 = checksums[3L],
    Chr22GFF3_MD5 = unname(tools::md5sum(gff))
  )
  write.dcf(as.data.frame(as.list(metadata)), file.path(output, "metadata.dcf"))
  print(comparisons, row.names = FALSE)
  if (any(comparisons$status != "equivalent") || !all(vcf$same_samples)) {
    stop("The readers differ from the source data. See the comparison reports.", call. = FALSE)
  }
}

if (sys.nframe() == 0L) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  script <- normalizePath(script, winslash = "/", mustWork = TRUE)
  source(file.path(dirname(script), "compare_standard_io.R"))
  source(file.path(dirname(script), "assess_standard_backends.R"))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 3L && args[1L] == "--worker" && args[2L] %in% names(real_io_readers)) {
    real_io_worker(args[2L], args[3L])
  } else if (length(args) == 2L) {
    assess_real_standard_io(args[1L], args[2L], script)
  } else {
    stop("Usage: Rscript tools/assess_real_standard_io.R DATA_DIR OUTPUT_DIR", call. = FALSE)
  }
}
