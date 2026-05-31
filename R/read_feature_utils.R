# Author: Rensc
# Date: 2026-05-28
# Version: 0.2.4
# Function: Shared parsers for feature annotation files
# Input: GFF/GTF attribute strings and feature tables
# Output: Standardized feature tables


make_progress_message <- function(verbose) {
  force(verbose)
  function(..., appendLF = TRUE) {
    if (isTRUE(verbose)) {
      message(sprintf(...), appendLF = appendLF)
    }
  }
}

make_stage_progress <- function(total, progress = TRUE, prefix = "[GeneTrackR]") {
  total <- as.integer(total)[1L]
  if (!isTRUE(progress) || is.na(total) || total <= 0L) {
    return(list(
      tick = function(value = 1L) invisible(NULL),
      close = function() invisible(NULL)
    ))
  }

  current <- 0L
  width <- 30L

  render <- function(value) {
    pct <- if (total > 0L) value / total else 1
    filled <- max(0L, min(width, round(width * pct)))
    bar <- paste0(strrep("=", filled), strrep(".", width - filled))
    sprintf(
      "%s Progress [%s] %s/%s (%3.0f%%)",
      prefix,
      bar,
      format(value, big.mark = ","),
      format(total, big.mark = ","),
      pct * 100
    )
  }

  list(
    tick = function(value = NULL) {
      if (is.null(value)) {
        current <<- min(total, current + 1L)
      } else {
        current <<- min(total, max(0L, as.integer(value)[1L]))
      }
      message(render(current))
      invisible(current)
    },
    close = function() invisible(NULL)
  )
}

with_quiet_datatable <- function(expr) {
  old <- options(datatable.verbose = FALSE)
  on.exit(options(old), add = TRUE)
  force(expr)
}

get_file_size_bytes <- function(file) {
  suppressWarnings(as.numeric(file.info(file)$size)[1L])
}

format_file_size <- function(size) {
  size <- suppressWarnings(as.numeric(size)[1L])
  if (!is.finite(size) || is.na(size)) {
    return("unknown size")
  }
  if (size < 0) {
    return("unknown size")
  }

  units <- c("B", "KB", "MB", "GB", "TB")
  unit_id <- 1L
  while (size >= 1024 && unit_id < length(units)) {
    size <- size / 1024
    unit_id <- unit_id + 1L
  }

  if (unit_id == 1L) {
    paste0(format(round(size), big.mark = ",", scientific = FALSE), " ", units[unit_id])
  } else {
    paste0(format(round(size, 2), nsmall = 2, big.mark = ",", scientific = FALSE), " ", units[unit_id])
  }
}

format_file_size_message <- function(prefix, file) {
  input_file <- normalizePath(file, winslash = "/", mustWork = FALSE)
  file_size <- get_file_size_bytes(file)
  if (is.finite(file_size) && !is.na(file_size)) {
    sprintf("[GeneTrackR] Reading %s file: %s (%s)", prefix, input_file, format_file_size(file_size))
  } else {
    sprintf("[GeneTrackR] Reading %s file: %s", prefix, input_file)
  }
}

read_gff_gtf <- function(file,
                         format = c("GFF", "GTF"),
                         feature_types = NULL,
                         verbose = TRUE,
                         progress = interactive() && isTRUE(verbose)) {
  format <- match.arg(format)
  verbose <- isTRUE(verbose)
  progress <- isTRUE(progress)
  old_dt_options <- options(datatable.verbose = FALSE)
  on.exit(options(old_dt_options), add = TRUE)
  stop_if_not(file.exists(file), paste0("File does not exist: ", file))

  progress_msg <- make_progress_message(verbose)
  stage <- make_stage_progress(total = 5L, progress = progress)
  on.exit(stage$close(), add = TRUE)

  input_file <- normalizePath(file, winslash = "/", mustWork = FALSE)
  progress_msg("%s", format_file_size_message(format, file))

  dt <- data.table::fread(
    file,
    header = FALSE,
    sep = "\t",
    data.table = TRUE,
    comment.char = "#",
    showProgress = FALSE
  )
  stage$tick()
  stop_if_not(ncol(dt) >= 9L, paste0("A ", format, " file must contain at least 9 columns."))
  data.table::setnames(dt, seq_len(9L), c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attribute"))
  progress_msg("[GeneTrackR] Loaded %s records with %s columns.", format(nrow(dt), big.mark = ","), ncol(dt))
  progress_msg("[GeneTrackR] Detected input format: %s", format)

  if (!is.null(feature_types)) {
    n_before <- nrow(dt)
    dt <- dt[as.character(dt[["type"]]) %in% as.character(feature_types)]
    progress_msg(
      "[GeneTrackR] Kept %s/%s records after feature type filtering.",
      format(nrow(dt), big.mark = ","),
      format(n_before, big.mark = ",")
    )
  }
  stage$tick()

  progress_msg("[GeneTrackR] Parsing %s attribute column.", format)
  attrs <- parse_feature_attributes_fast(dt[["attribute"]], format = format)
  stage$tick()

  progress_msg("[GeneTrackR] Building standardized feature table.")
  dt[, "parent_id" := pick_first_nonempty(attrs$Parent, NA_character_)]
  dt[, "feature_id" := pick_first_nonempty(attrs$ID, attrs$transcript_id, attrs$gene_id, attrs$Name, paste0(format, "_", seq_len(.N)))]
  dt[, "name" := pick_first_nonempty(attrs$Name, attrs$gene_name, attrs$gene_id, attrs$transcript_id, dt[["feature_id"]])]
  dt[, "level" := infer_feature_level(type)]
  dt[, "gene_id" := pick_first_nonempty(attrs$gene_id, data.table::fifelse(level == "gene", feature_id, NA_character_), NA_character_)]
  dt[, "transcript_id" := pick_first_nonempty(attrs$transcript_id, data.table::fifelse(level == "transcript", feature_id, NA_character_), data.table::fifelse(level %in% c("exon", "subfeature"), parent_id, NA_character_), NA_character_)]
  dt[level == "transcript" & (is.na(gene_id) | gene_id == ""), "gene_id" := parent_id]
  dt[level %in% c("exon", "subfeature") & (is.na(gene_id) | gene_id == ""), "gene_id" := attrs$gene_id]
  dt[, "gene_type" := pick_first_nonempty(attrs$gene_biotype, attrs$gene_type, attrs$transcript_biotype, attrs$transcript_type, attrs$biotype, NA_character_)]
  dt[, "exon_number" := suppressWarnings(as.integer(pick_first_nonempty(attrs$exon_number, NA_character_)))]
  dt[, "score" := suppressWarnings(as.numeric(dt[["score"]]))]
  dt[is.na(score), "score" := NA_real_]

  out <- dt[, .(
    feature_id = as.character(feature_id),
    name = as.character(name),
    chrom = as.character(chrom),
    start = as.integer(start),
    end = as.integer(end),
    type = as.character(type),
    level = as.character(level),
    score = as.numeric(score),
    strand = as.character(strand),
    source = as.character(source),
    gene_id = as.character(gene_id),
    transcript_id = as.character(transcript_id),
    parent_id = as.character(parent_id),
    gene_type = as.character(gene_type),
    exon_number = as.integer(exon_number),
    phase = as.character(phase),
    attribute = as.character(attribute)
  )]
  stage$tick()

  progress_msg("[GeneTrackR] Building Feature object and derived gene/transcript/exon tables.")
  hierarchy <- derive_feature_hierarchy_fast(out)
  obj <- FeatureTrack(
    out,
    genes = hierarchy$genes,
    transcripts = hierarchy$transcripts,
    exons = hierarchy$exons,
    meta = list(
      source_file = input_file,
      format = format,
      coordinate_input = "1-based closed",
      coordinate_internal = "1-based closed"
    )
  )
  stage$tick()
  progress_msg(
    "[GeneTrackR] Finished. records: %s; genes: %s; transcripts: %s; exons: %s.",
    format(nrow(obj$data), big.mark = ","),
    format(if (!is.null(obj$genes)) nrow(obj$genes) else 0L, big.mark = ","),
    format(if (!is.null(obj$transcripts)) nrow(obj$transcripts) else 0L, big.mark = ","),
    format(if (!is.null(obj$exons)) nrow(obj$exons) else 0L, big.mark = ",")
  )
  obj
}

parse_feature_attributes_fast <- function(x, format = c("GFF", "GTF"), keys = NULL) {
  format <- match.arg(format)
  if (is.null(keys)) {
    keys <- c(
      "ID", "Name", "Parent", "gene_id", "transcript_id", "gene_name",
      "gene_biotype", "gene_type", "transcript_biotype", "transcript_type",
      "biotype", "exon_number", "protein_id", "ccds_id"
    )
  }
  out <- stats::setNames(vector("list", length(keys)), keys)
  x <- as.character(x)
  x[is.na(x)] <- ""
  for (key in keys) {
    out[[key]] <- extract_feature_attribute(x, key = key, format = format)
  }
  out
}

parse_feature_attributes <- parse_feature_attributes_fast

extract_feature_attribute <- function(x, key, format = c("GFF", "GTF")) {
  format <- match.arg(format)
  x <- as.character(x)
  x[is.na(x)] <- ""
  key_pat <- gsub("([\\.\\+\\*\\?\\[\\^\\]\\$\\(\\)\\{\\}\\=\\!\\<\\>\\|\\:\\-])", "\\\\\\1", key, perl = TRUE)
  if (format == "GFF") {
    pattern <- paste0("(?:^|;)\\s*", key_pat, "=([^;]+)")
    replace_pattern <- paste0("^.*?", key_pat, "=([^;]+).*$")
  } else {
    pattern <- paste0("(?:^|;)\\s*", key_pat, "\\s+\\\"?([^\\\";]+)\\\"?")
    replace_pattern <- paste0("^.*?", key_pat, "\\s+\\\"?([^\\\";]+)\\\"?.*$")
  }
  m <- regexpr(pattern, x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  hit <- which(m > 0L)
  if (length(hit) > 0L) {
    matched <- regmatches(x, m)[seq_along(hit)]
    out[hit] <- trimws(sub(replace_pattern, "\\1", matched, perl = TRUE))
    out[out == ""] <- NA_character_
  }
  out
}

pick_first_nonempty <- function(...) {
  vals <- list(...)
  n <- max(vapply(vals, length, integer(1L)))
  out <- rep(NA_character_, n)
  for (v in vals) {
    vv <- rep(as.character(v), length.out = n)
    idx <- is.na(out) | !nzchar(out)
    out[idx & !is.na(vv) & nzchar(vv)] <- vv[idx & !is.na(vv) & nzchar(vv)]
  }
  out
}


derive_feature_hierarchy_fast <- function(dt) {
  x <- data.table::copy(data.table::as.data.table(dt))
  data.table::setorderv(x, c("chrom", "start", "end", "feature_id"))

  gene_rows <- x[tolower(x[["level"]]) == "gene" | tolower(x[["type"]]) == "gene"]
  tx_rows <- x[tolower(x[["level"]]) == "transcript" | is_transcript_feature_type(x[["type"]])]
  exon_rows <- x[tolower(x[["level"]]) == "exon" | tolower(x[["type"]]) == "exon"]
  cds_rows <- x[tolower(x[["type"]]) == "cds"]

  if (nrow(tx_rows) > 0L) {
    tx <- data.table::copy(tx_rows)
    tx[is.na(transcript_id) | transcript_id == "", "transcript_id" := feature_id]
    tx[is.na(gene_id) | gene_id == "", "gene_id" := parent_id]
    tx[is.na(gene_id) | gene_id == "", "gene_id" := transcript_id]
    transcripts <- tx[, .(
      transcript_id = as.character(transcript_id),
      gene_id = as.character(gene_id),
      chrom = as.character(chrom),
      strand = as.character(strand),
      tx_start = as.integer(start),
      tx_end = as.integer(end),
      gene_type = as.character(gene_type),
      score = as.numeric(score)
    )]
  } else if (nrow(exon_rows) > 0L && any(!is.na(exon_rows[["transcript_id"]]))) {
    transcripts <- exon_rows[!is.na(transcript_id), .(
      gene_id = pick_first_nonempty(gene_id[1L], parent_id[1L], transcript_id[1L]),
      chrom = chrom[1L],
      strand = strand[1L],
      tx_start = as.integer(min(start, na.rm = TRUE)),
      tx_end = as.integer(max(end, na.rm = TRUE)),
      gene_type = NA_character_,
      score = NA_real_
    ), by = transcript_id]
  } else {
    transcripts <- data.table::data.table()
  }

  if (nrow(transcripts) > 0L) {
    cds_summary <- data.table::data.table(transcript_id = character(), cds_start = integer(), cds_end = integer())
    if (nrow(cds_rows) > 0L) {
      cds_rows <- data.table::copy(cds_rows)
      cds_rows[is.na(transcript_id) | transcript_id == "", "transcript_id" := parent_id]
      cds_summary <- cds_rows[!is.na(transcript_id), .(
        cds_start = as.integer(min(start, na.rm = TRUE)),
        cds_end = as.integer(max(end, na.rm = TRUE))
      ), by = transcript_id]
    }
    transcripts <- merge(transcripts, cds_summary, by = "transcript_id", all.x = TRUE)
    transcripts[is.na(cds_start) | is.na(cds_end), `:=`(cds_start = as.integer(tx_start), cds_end = as.integer(tx_start - 1L))]
    transcripts[, "gene_type" := data.table::fifelse(as.integer(cds_start) <= as.integer(cds_end), "coding", "non-coding")]
    transcripts[, "exon_count" := NA_integer_]
    transcripts[, "cds_start_stat" := data.table::fifelse(as.character(gene_type) == "coding", "cmpl", "none")]
    transcripts[, "cds_end_stat" := data.table::fifelse(as.character(gene_type) == "coding", "cmpl", "none")]
    transcripts[, "exon_frames" := NA_character_]
  }

  if (nrow(exon_rows) > 0L) {
    exons <- data.table::copy(exon_rows)
    exons[is.na(transcript_id) | transcript_id == "", "transcript_id" := parent_id]
    if (nrow(transcripts) > 0L) {
      tx_map <- transcripts[, .(transcript_id, gene_id)]
      exons <- merge(exons, tx_map, by = "transcript_id", all.x = TRUE, suffixes = c("", ".tx"))
      exons[is.na(gene_id) | gene_id == "", "gene_id" := gene_id.tx]
      exons[, "gene_id.tx" := NULL]
    }
    exons[is.na(gene_id) | gene_id == "", "gene_id" := parent_id]
    data.table::setorderv(exons, c("transcript_id", "start", "end"))
    exons[is.na(exon_number), "exon_number" := seq_len(.N), by = transcript_id]
    exons <- exons[, .(
      transcript_id = as.character(transcript_id),
      gene_id = as.character(gene_id),
      chrom = as.character(chrom),
      strand = as.character(strand),
      exon_number = as.integer(exon_number),
      exon_start = as.integer(start),
      exon_end = as.integer(end),
      exon_frame = -1L
    )]
  } else {
    exons <- data.table::data.table()
  }

  if (nrow(transcripts) > 0L && nrow(exons) > 0L) {
    exon_counts <- exons[, .(exon_count = as.integer(.N)), by = transcript_id]
    transcripts[, "exon_count" := NULL]
    transcripts <- merge(transcripts, exon_counts, by = "transcript_id", all.x = TRUE)
    transcripts[is.na(exon_count), "exon_count" := 0L]
  }

  frame_info <- fill_cds_status_and_exon_frames(transcripts, exons, cds_rows)
  transcripts <- frame_info$transcripts
  exons <- frame_info$exons

  if (nrow(transcripts) > 0L) {
    genes <- build_gene_table(transcripts)
  } else if (nrow(gene_rows) > 0L) {
    genes <- gene_rows[, .(
      gene_id = pick_first_nonempty(gene_id, feature_id, name),
      chrom = chrom,
      strand = strand,
      gene_start = as.integer(start),
      gene_end = as.integer(end),
      n_transcripts = 0L,
      gene_type = data.table::fifelse(!is.na(gene_type) & gene_type != "", gene_type, "unknown")
    )]
  } else {
    genes <- data.table::data.table()
  }

  list(genes = genes, transcripts = transcripts, exons = exons)
}
