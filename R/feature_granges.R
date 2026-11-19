# Author: Rensc
# Date: 2026-11-19
# Version: dev002
# Function: Adapt GRanges and Feature objects at the interval-feature level
# Input: GRanges or Feature-compatible annotation objects
# Output: Feature or GRanges objects with aligned interval metadata

granges_feature_column <- function(value, name) {
  stop_if_not(
    is.null(dim(value)),
    sprintf("GRanges metadata column '%s' must be one-dimensional.", name)
  )
  if (inherits(value, "Rle")) {
    value <- S4Vectors::decode(value)
  } else if (inherits(value, "List")) {
    value <- as.list(value)
  } else {
    stop_if_not(
      !isS4(value) && (is.atomic(value) || is.list(value)),
      sprintf("GRanges metadata column '%s' has an unsupported class.", name)
    )
  }

  character_fields <- c(
    "feature_id", "name", "type", "level", "source", "gene_id",
    "transcript_id", "parent_id", "gene_type", "phase", "attribute"
  )
  if (name %in% c(character_fields, "score", "exon_number")) {
    stop_if_not(
      is.atomic(value),
      sprintf("GRanges metadata column '%s' must be an atomic vector for Feature annotations.", name)
    )
  }
  if (name %in% c("score", "exon_number")) {
    if (is.factor(value)) value <- as.character(value)
    if (is.character(value)) value <- utils::type.convert(value, as.is = TRUE)
    stop_if_not(
      is.numeric(value) || all(is.na(value)),
      sprintf("GRanges metadata column '%s' must be numeric or entirely missing.", name)
    )
  }
  if (identical(name, "exon_number") && !all(is.na(value))) {
    stop_if_not(
      all(is.na(value) | (
        is.finite(value) & abs(value) <= .Machine$integer.max & value == trunc(value)
      )),
      "GRanges metadata column 'exon_number' must contain integers or missing values."
    )
  }
  data.table::copy(value)
}

#' @rdname as_feature
#' @export
as_feature.GRanges <- function(object) {
  methods::validObject(object)
  start <- as.integer(IRanges::start(object))
  end <- as.integer(IRanges::end(object))
  stop_if_not(
    all(start >= 1L),
    "GRanges input must use positive, 1-based start coordinates."
  )
  stop_if_not(
    all(end >= start),
    "GRanges input must contain positive-width intervals; zero-width ranges cannot represent Feature intervals."
  )
  metadata <- S4Vectors::mcols(object)
  metadata_names <- names(metadata)
  stop_if_not(
    !anyNA(metadata_names) && all(nzchar(metadata_names)) && !anyDuplicated(metadata_names),
    "GRanges metadata column names must be unique and nonempty."
  )
  structural <- intersect(metadata_names, c(
    "chrom", "start", "end", "strand", ".N", ".I", ".SD", ".BY", ".GRP"
  ))
  stop_if_not(
    length(structural) == 0L,
    sprintf("GRanges metadata conflicts with Feature structural columns: %s. Rename these metadata columns before conversion.", paste(structural, collapse = ", "))
  )

  data <- data.table::data.table(
    chrom = as.character(GenomicRanges::seqnames(object)),
    start = start,
    end = end,
    strand = as.character(GenomicRanges::strand(object))
  )
  for (name in metadata_names) {
    value <- granges_feature_column(metadata[[name]], name)
    data.table::set(data, j = name, value = list(value))
  }

  if (!is.null(names(object)) && !"name" %in% names(data)) {
    data.table::set(data, j = "name", value = list(names(object)))
  }
  row_column <- utils::tail(make.unique(c(names(data), ".gtr_granges_row_id"), sep = "_"), 1L)
  data.table::set(data, j = row_column, value = list(seq_len(length(object))))
  Feature(
    data = data,
    meta = list(
      format = "GRanges",
      coordinate_internal = "1-based closed",
      granges_context = list(
        origin = "GRanges",
        seqinfo = GenomicRanges::seqinfo(object),
        metadata = S4Vectors::metadata(object),
        source_mcols = metadata,
        range_names = names(object),
        row_column = row_column
      )
    )
  )
}

restore_granges_feature_column <- function(source, rows, value, standard) {
  selected <- source[rows]
  if (standard && !isS4(selected)) {
    names(value) <- names(selected)
    return(value)
  }
  if (inherits(selected, "Vector")) {
    annotations <- S4Vectors::mcols(selected)
    metadata <- S4Vectors::metadata(selected)
  }
  selected_names <- names(selected)
  selected[] <- value
  if (inherits(selected, "Vector")) {
    S4Vectors::mcols(selected) <- annotations
    S4Vectors::metadata(selected) <- metadata
  }
  names(selected) <- selected_names
  selected
}

feature_granges_metadata <- function(data, context, rows) {
  standard <- c(
    "feature_id", "name", "type", "gene_id", "transcript_id", "level", "score",
    "source", "parent_id", "gene_type", "exon_number", "phase", "attribute"
  )
  columns <- unique(c(
    intersect(standard[seq_len(5L)], names(data)),
    setdiff(names(data), c("chrom", "start", "end", "strand", context$row_column))
  ))
  reserved <- intersect(columns, c(
    "seqnames", "ranges", "strand", "seqlevels", "seqlengths",
    "isCircular", "start", "end", "width", "element"
  ))
  stop_if_not(length(reserved) == 0L, sprintf(
    "Feature metadata conflicts with reserved GRanges column names: %s. Rename these columns before conversion.",
    paste(reserved, collapse = ", ")
  ))
  metadata <- S4Vectors::make_zero_col_DFrame(nrow(data))
  for (name in columns) {
    value <- data.table::copy(data[[name]])
    if (name %in% names(context$source_mcols)) {
      value <- restore_granges_feature_column(
        context$source_mcols[[name]], rows, value, name %in% standard
      )
    }
    metadata[[name]] <- value
  }
  if (!is.null(context)) {
    source <- context$source_mcols
    S4Vectors::metadata(metadata) <- S4Vectors::metadata(source)
    annotations <- S4Vectors::mcols(source)
    if (!is.null(annotations)) {
      annotations <- annotations[match(names(metadata), names(source)), , drop = FALSE]
      rownames(annotations) <- names(metadata)
      S4Vectors::mcols(metadata) <- annotations
    }
  }
  metadata
}

feature_as_granges <- function(object) {
  data <- as_feature_table(object)
  stop_if_not(
    !anyDuplicated(names(data)),
    "Feature column names must be unique before conversion to GRanges."
  )
  context <- object$meta$granges_context
  if (!is.list(context) || !identical(context$origin, "GRanges")) {
    context <- NULL
  }
  rows <- NULL
  if (!is.null(context)) {
    required <- c(names(context$source_mcols), context$row_column)
    stop_if_not(
      all(required %in% names(data)),
      "GRanges names or metadata were lost during a derived Feature operation. Subset the GRanges before conversion when rebuilding gene models."
    )
    rows <- data[[context$row_column]]
    stop_if_not(
      is.numeric(rows) && !anyNA(rows) && all(
        is.finite(rows) & rows == trunc(rows) & rows >= 1L & rows <= nrow(context$source_mcols)
      ),
      "GRanges source-row indices must remain valid integers after Feature operations."
    )
  }
  ranges <- GenomicRanges::GRanges(
    seqnames = data$chrom,
    ranges = IRanges::IRanges(start = data$start, end = data$end),
    strand = data$strand,
    seqinfo = context$seqinfo
  )
  S4Vectors::mcols(ranges) <- feature_granges_metadata(data, context, rows)
  if (!is.null(context)) {
    S4Vectors::metadata(ranges) <- context$metadata
    if (!is.null(context$range_names)) names(ranges) <- context$range_names[rows]
  }
  methods::validObject(ranges)
  ranges
}
