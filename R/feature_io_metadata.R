# Author: Rensc
# Date: 2026-11-20
# Version: dev001
# Function: Preserve and validate format-specific annotation metadata.
# Input: GFF3 parent tokens or flat Feature BED12 metadata.
# Output: Structured metadata and validated standard-format output fields.

gff_parent_ids <- function(values) {
    lapply(as.character(values), function(value) {
        if (is.na(value) || !nzchar(value) || value == ".") return(character())
        parents <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
        stop_if_not(all(nzchar(parents)), "GFF3 Parent contains empty identifiers.")
        parents
    })
}

gff_parent_attribute <- function(data) {
    primary <- as.character(data$parent_id)
    if (!"parent_ids" %in% names(data)) return(primary)
    parents <- data$parent_ids
    if (inherits(parents, "List")) parents <- as.list(parents)
    stop_if_not(is.list(parents) && length(parents) == nrow(data),
        "`parent_ids` must be a row-aligned character-list column.")
    vapply(seq_len(nrow(data)), function(i) {
        values <- parents[[i]]
        if (is.null(values)) return(primary[i])
        stop_if_not(is.character(values) && !anyNA(values) &&
            all(nzchar(values)) && !any(grepl("[,;\t\r\n]", values)),
            "`parent_ids` must contain nonempty raw GFF3 identifier tokens.")
        first <- if (length(values)) values[1L] else NA_character_
        current <- if (is.na(primary[i]) || !nzchar(primary[i])) NA_character_ else primary[i]
        stop_if_not(identical(unname(current), unname(first)),
            "`parent_id` must match the first `parent_ids` entry; update both fields consistently.")
        if (!length(values)) NA_character_ else paste(values, collapse = ",")
    }, character(1))
}

guard_gff_parent_rebuild <- function(object, transcripts) {
    data <- object$data
    if (!"parent_ids" %in% names(data)) return(invisible(NULL))
    ids <- unique(as.character(c(transcripts$transcript_id, transcripts$gene_id)))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    selected <- as.character(data$transcript_id) %in% ids |
        as.character(data$gene_id) %in% ids | as.character(data$feature_id) %in% ids
    linked <- vapply(data$parent_ids, function(parents) any(parents %in% ids), logical(1))
    stop_if_not(!any((selected | linked) & lengths(data$parent_ids) > 1L),
        paste("Rebuilding the primary gene model would discard complete GFF3 parent links.",
            "Subset feature-level GRanges before converting back with as_feature()."))
    invisible(NULL)
}

retrieve_primary_feature <- function(object, ...) {
    # Internal calculations intentionally consume the established primary model.
    # They do not export a replacement for the original annotation metadata.
    projected <- object
    if (is_gene_model_feature(object) && "parent_ids" %in% names(object$data)) {
        projected$data <- data.table::copy(object$data)
        projected$data[["parent_ids"]] <- NULL
        projected$meta$granges_context <- NULL
    }
    retrieve_feature(projected, ...)
}

bed12_integer <- function(value, field, minimum = 0L) {
    if (is.character(value)) {
        stop_if_not(!anyNA(value) && all(grepl("^[0-9]+$", value)),
            paste0("BED12 ", field, " must contain exact integers."))
        value <- as.numeric(value)
    }
    stop_if_not(is.numeric(value) && !is.factor(value) && !anyNA(value) &&
        all(is.finite(value) & value >= minimum &
            value <= .Machine$integer.max & value == floor(value)),
        paste0("BED12 ", field, " contains invalid or out-of-range integers."))
    as.integer(value)
}

bed12_block_list <- function(value, field, minimum) {
    lapply(as.character(value), function(text) {
        stop_if_not(!is.na(text) && grepl("^[0-9]+(,[0-9]+)*,?$", text),
            paste0("BED12 ", field, " must be a comma-separated integer list."))
        bed12_integer(strsplit(sub(",$", "", text), ",", fixed = TRUE)[[1L]],
            field, minimum)
    })
}

validate_bed12_metadata <- function(data, metadata) {
    rows <- nrow(data)
    required <- c("thick_start", "thick_end", "item_rgb", "block_count",
        "block_sizes", "block_starts")
    stop_if_not(all(required %in% names(metadata)), "BED12 metadata is incomplete: missing fields.")
    stop_if_not(all(vapply(metadata[required], length, integer(1)) == rows),
        "BED12 metadata must be row-aligned.")
    start <- as.numeric(bed12_integer(data$start, "start", 1L))
    end <- as.numeric(bed12_integer(data$end, "end", 1L))
    stop_if_not(all(end >= start), "BED12 intervals must have positive width.")
    metadata$thick_start <- bed12_integer(metadata$thick_start, "thick_start", 1L)
    metadata$thick_end <- bed12_integer(metadata$thick_end, "thick_end")
    metadata$block_count <- bed12_integer(metadata$block_count, "block_count", 1L)
    thick_start <- as.numeric(metadata$thick_start)
    thick_end <- as.numeric(metadata$thick_end)
    stop_if_not(all(thick_start >= start & thick_start <= end + 1 &
        thick_end >= start - 1 & thick_end <= end & thick_end >= thick_start - 1),
        "BED12 thick intervals must be contained, or explicitly empty.")
    for (field in c("block_sizes", "block_starts")) {
        if (inherits(metadata[[field]], "List")) metadata[[field]] <- as.list(metadata[[field]])
        stop_if_not(is.list(metadata[[field]]), paste0("BED12 ", field, " must be a list column."))
        metadata[[field]] <- lapply(metadata[[field]], bed12_integer,
            field = field, minimum = if (field == "block_sizes") 1L else 0L)
    }
    validate_bed12_blocks(metadata, end - start + 1)
    metadata$item_rgb <- as.character(metadata$item_rgb)
    valid_rgb <- vapply(metadata$item_rgb, function(value) {
        if (is.na(value)) return(FALSE)
        if (value == "0") return(TRUE)
        if (!grepl("^[0-9]+,[0-9]+,[0-9]+$", value)) return(FALSE)
        all(as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]]) <= 255)
    }, logical(1))
    stop_if_not(all(valid_rgb), "BED12 item_rgb must be 0 or three RGB bytes.")
    metadata
}

validate_bed12_blocks <- function(metadata, widths) {
    for (i in seq_along(widths)) {
        sizes <- as.numeric(metadata$block_sizes[[i]])
        offsets <- as.numeric(metadata$block_starts[[i]])
        count <- metadata$block_count[i]
        stop_if_not(length(sizes) == count && length(offsets) == count,
            sprintf("BED12 block counts disagree at row %d.", i))
        stop_if_not(offsets[1L] == 0 && all(offsets + sizes <= widths[i]) &&
            utils::tail(offsets + sizes, 1L) == widths[i] &&
            (count == 1L || all(offsets[-1L] >= utils::head(offsets + sizes, -1L))),
            sprintf("BED12 block geometry is invalid at row %d.", i))
    }
    invisible(NULL)
}

validate_bed12_coordinates <- function(data, coordinate) {
    first <- bed12_integer(data$bed_start, "chromStart",
        if (coordinate == "bed") 0L else 1L)
    last <- bed12_integer(data$bed_end, "chromEnd", 1L)
    if (coordinate == "bed") {
        stop_if_not(all(first < .Machine$integer.max & last > first),
            "BED12 half-open intervals must have positive representable width.")
    } else {
        stop_if_not(all(last >= first), "BED12 closed intervals must have positive width.")
    }
    invisible(NULL)
}

read_bed12_metadata <- function(data, coordinate) {
    metadata <- list(
        thick_start = bed12_integer(data$V7, "thick_start"),
        thick_end = bed12_integer(data$V8, "thick_end"),
        item_rgb = as.character(data$V9),
        block_count = bed12_integer(data$V10, "block_count", 1L),
        block_sizes = bed12_block_list(data$V11, "block_sizes", 1L),
        block_starts = bed12_block_list(data$V12, "block_starts", 0L)
    )
    if (coordinate == "bed") {
        metadata$thick_start <- as.numeric(metadata$thick_start) + 1
    }
    validate_bed12_metadata(data, metadata)
}

has_bed12_metadata <- function(object) {
    !inherits(object, "GenePred") && (isTRUE(object$meta$bed12) ||
        any(c("block_count", "block_sizes", "block_starts") %in% names(object$data)))
}

bed12_feature_table <- function(object, coordinate, sort_output, chrom_order) {
    data <- data.table::copy(as_feature_table(object))
    fields <- c("thick_start", "thick_end", "item_rgb", "block_count",
        "block_sizes", "block_starts")
    stop_if_not(all(fields %in% names(data)), "BED12 metadata is incomplete: missing fields.")
    metadata <- validate_bed12_metadata(data, as.list(data[, fields, with = FALSE]))
    out <- data.table::data.table(
        chrom = as.character(data$chrom),
        chromStart = if (coordinate == "ucsc") as.integer(data$start) - 1L else as.integer(data$start),
        chromEnd = as.integer(data$end),
        name = ifelse(is.na(data$name) | !nzchar(data$name), ".", as.character(data$name)),
        score = normalize_bed_score(data$score), strand = normalize_bed_strand(data$strand),
        thickStart = if (coordinate == "ucsc") metadata$thick_start - 1L else metadata$thick_start,
        thickEnd = metadata$thick_end, itemRgb = metadata$item_rgb,
        blockCount = metadata$block_count,
        blockSizes = vapply(metadata$block_sizes, paste, character(1), collapse = ","),
        blockStarts = vapply(metadata$block_starts, paste, character(1), collapse = ",")
    )
    stop_if_not(!anyNA(out$chrom) && all(nzchar(out$chrom)) &&
        !any(grepl("[\t\r\n]", c(out$chrom, out$name))),
        "BED12 chromosome and name fields contain invalid text.")
    if (isTRUE(sort_output)) out <- sort_bed_output_table(out, chrom_order)
    out
}
