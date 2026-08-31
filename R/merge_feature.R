# Author: Rensc
# Date: 2026-09-01
# Version: dev005
# Function: Merge Feature/GenePred-compatible annotation objects
# Input: Feature-compatible annotation objects
# Output: Merged Feature/GenePred-compatible object


is_feature_merge_object <- function(object) {
  inherits(object, "Feature") ||
    inherits(object, "FeatureTrack") ||
    inherits(object, "GenePred")
}


normalize_feature_merge_inputs <- function(objects) {
  if (
    length(objects) == 1L &&
      is.list(objects[[1L]]) &&
      !is_feature_merge_object(objects[[1L]])
  ) {
    objects <- objects[[1L]]
  }
  objects
}


nonempty_merge_id <- function(x) {
  x <- as.character(x)
  !is.na(x) & nzchar(x)
}


collect_merge_id_sources <- function(tables, id_col) {
  pairs <- lapply(tables, function(dt) {
    if (
      is.null(dt) ||
        nrow(dt) == 0L ||
        !id_col %in% names(dt) ||
        !".merge_source_index" %in% names(dt)
    ) {
      return(NULL)
    }

    value <- as.character(dt[[id_col]])
    keep <- nonempty_merge_id(value)
    if (!any(keep)) {
      return(NULL)
    }

    data.table::data.table(
      merge_id = value[keep],
      merge_source_index = as.integer(dt[[".merge_source_index"]][keep])
    )
  })

  pairs <- data.table::rbindlist(pairs, use.names = TRUE, fill = TRUE)
  if (nrow(pairs) == 0L) {
    return(data.table::data.table(
      merge_id = character(),
      first_source_index = integer(),
      source_count = integer()
    ))
  }

  pairs <- unique(pairs)
  pairs[, .(
    first_source_index = min(merge_source_index),
    source_count = data.table::uniqueN(merge_source_index)
  ), by = merge_id]
}


make_merge_id_source_table <- function(dt, id_col) {
  if (
    is.null(dt) ||
      nrow(dt) == 0L ||
      !id_col %in% names(dt) ||
      !".merge_source_index" %in% names(dt)
  ) {
    return(data.table::data.table(
      merge_id = character(),
      .merge_source_index = integer()
    ))
  }

  data.table::data.table(
    merge_id = as.character(dt[[id_col]]),
    .merge_source_index = as.integer(dt[[".merge_source_index"]])
  )
}


format_merge_duplicate_warning <- function(duplicate_ids, conflict) {
  parts <- vapply(names(duplicate_ids), function(id_type) {
    ids <- duplicate_ids[[id_type]]
    if (length(ids) == 0L) {
      return(NA_character_)
    }
    examples <- paste(utils::head(ids, 5L), collapse = ", ")
    if (length(ids) > 5L) {
      examples <- paste0(examples, ", ...")
    }
    sprintf("%s=%s [%s]", id_type, format(length(ids), big.mark = ","), examples)
  }, character(1L))
  parts <- parts[!is.na(parts)]

  sprintf(
    paste0(
      "Duplicated annotation identifiers were detected across input objects: %s. ",
      "Applying `conflict = \"%s\"`. With `conflict = \"deduplicate\"`, input ",
      "order determines which record is retained."
    ),
    paste(parts, collapse = "; "),
    conflict
  )
}


filter_merge_id_owner <- function(dt, id_col, ownership) {
  if (
    is.null(dt) ||
      nrow(dt) == 0L ||
      !id_col %in% names(dt) ||
      nrow(ownership) == 0L
  ) {
    return(dt)
  }

  owner <- stats::setNames(
    ownership[["first_source_index"]],
    ownership[["merge_id"]]
  )
  value <- as.character(dt[[id_col]])
  expected_source <- unname(owner[value])
  keep <- is.na(expected_source) |
    as.integer(dt[[".merge_source_index"]]) == as.integer(expected_source)
  dt[keep]
}


filter_merge_gene_rows <- function(dt, ownership) {
  if (
    is.null(dt) ||
      nrow(dt) == 0L ||
      !"gene_id" %in% names(dt) ||
      nrow(ownership) == 0L
  ) {
    return(dt)
  }

  is_gene_row <- rep(FALSE, nrow(dt))
  if ("level" %in% names(dt)) {
    is_gene_row <- is_gene_row |
      tolower(as.character(dt[["level"]])) == "gene"
  }
  if ("type" %in% names(dt)) {
    is_gene_row <- is_gene_row |
      tolower(as.character(dt[["type"]])) == "gene"
  }
  if (!any(is_gene_row)) {
    return(dt)
  }

  owner <- stats::setNames(
    ownership[["first_source_index"]],
    ownership[["merge_id"]]
  )
  value <- as.character(dt[["gene_id"]])
  expected_source <- unname(owner[value])
  remove <- is_gene_row &
    !is.na(expected_source) &
    as.integer(dt[[".merge_source_index"]]) != as.integer(expected_source)
  dt[!remove]
}


update_merged_gene_feature_rows <- function(data, genes) {
  if (
    is.null(data) ||
      nrow(data) == 0L ||
      is.null(genes) ||
      nrow(genes) == 0L ||
      !"gene_id" %in% names(data)
  ) {
    return(data)
  }

  is_gene_row <- rep(FALSE, nrow(data))
  if ("level" %in% names(data)) {
    is_gene_row <- is_gene_row |
      tolower(as.character(data[["level"]])) == "gene"
  }
  if ("type" %in% names(data)) {
    is_gene_row <- is_gene_row |
      tolower(as.character(data[["type"]])) == "gene"
  }
  row_index <- which(is_gene_row)
  if (length(row_index) == 0L) {
    return(data)
  }

  gene_keys <- as.character(genes[["gene_id"]])
  data_keys <- as.character(data[["gene_id"]][row_index])
  gene_index <- match(data_keys, gene_keys)
  matched <- !is.na(gene_index)
  if (!any(matched)) {
    return(data)
  }

  assignments <- list(
    chrom = "chrom",
    start = "gene_start",
    end = "gene_end",
    strand = "strand",
    gene_type = "gene_type"
  )
  for (data_col in names(assignments)) {
    gene_col <- assignments[[data_col]]
    if (!data_col %in% names(data) || !gene_col %in% names(genes)) {
      next
    }
    data.table::set(
      data,
      i = row_index[matched],
      j = data_col,
      value = genes[[gene_col]][gene_index[matched]]
    )
  }
  data
}


deduplicate_feature_merge <- function(data,
                                      genes,
                                      transcripts,
                                      exons,
                                      gene_ownership,
                                      transcript_ownership,
                                      feature_ownership) {
  # A duplicated transcript is treated as one complete model. All exon and
  # subfeature records from later copies are removed together.
  transcripts <- filter_merge_id_owner(
    transcripts,
    "transcript_id",
    transcript_ownership
  )
  exons <- filter_merge_id_owner(
    exons,
    "transcript_id",
    transcript_ownership
  )
  data <- filter_merge_id_owner(
    data,
    "transcript_id",
    transcript_ownership
  )

  # Gene IDs may legitimately connect different, non-duplicated transcripts.
  # Only duplicate gene-level rows are removed; gene ranges are rebuilt below.
  genes <- filter_merge_id_owner(genes, "gene_id", gene_ownership)
  data <- filter_merge_gene_rows(data, gene_ownership)

  # Flat feature IDs identify individual annotation records.
  data <- filter_merge_id_owner(data, "feature_id", feature_ownership)

  if (!is.null(transcripts) && nrow(transcripts) > 0L) {
    rebuilt_genes <- build_gene_table(transcripts)

    if ("track_source" %in% names(transcripts)) {
      source_map <- transcripts[, .(
        track_source = paste(unique(as.character(track_source)), collapse = ",")
      ), by = gene_id]
      rebuilt_genes <- merge(
        rebuilt_genes,
        source_map,
        by = "gene_id",
        all.x = TRUE,
        sort = FALSE
      )
    }

    if (!is.null(genes) && nrow(genes) > 0L) {
      gene_only <- genes[
        !as.character(gene_id) %in% as.character(rebuilt_genes[["gene_id"]])
      ]
      genes <- data.table::rbindlist(
        list(rebuilt_genes, gene_only),
        use.names = TRUE,
        fill = TRUE
      )
    } else {
      genes <- rebuilt_genes
    }
    data <- update_merged_gene_feature_rows(data, genes)
  }

  list(
    data = data,
    genes = genes,
    transcripts = transcripts,
    exons = exons
  )
}


make_feature_merge_rename_map <- function(id_sources,
                                          rename_prefix,
                                          all_existing_ids) {
  duplicated_sources <- id_sources[source_count > 1L]
  if (nrow(duplicated_sources) == 0L) {
    return(data.table::data.table(
      merge_id = character(),
      merge_source_index = integer(),
      renamed_id = character()
    ))
  }

  source_pairs <- duplicated_sources[, {
    source_index <- seq.int(first_source_index + 1L, length(rename_prefix))
    list(merge_source_index = source_index)
  }, by = .(merge_id, first_source_index)]

  # Remove source/ID combinations that do not actually occur. The initial
  # sequence above is only a compact candidate generator.
  source_pairs <- source_pairs[
    paste(merge_id, merge_source_index, sep = "\r") %in%
      paste(
        attr(id_sources, "source_pairs")[["merge_id"]],
        attr(id_sources, "source_pairs")[["merge_source_index"]],
        sep = "\r"
      )
  ]

  used_ids <- unique(as.character(all_existing_ids))
  used_ids <- used_ids[nonempty_merge_id(used_ids)]
  renamed_id <- character(nrow(source_pairs))

  for (i in seq_len(nrow(source_pairs))) {
    source_index <- source_pairs[["merge_source_index"]][i]
    candidate <- paste0(
      as.character(rename_prefix[source_index]),
      as.character(source_pairs[["merge_id"]][i])
    )
    unique_candidate <- candidate
    suffix <- 2L
    while (unique_candidate %in% used_ids) {
      unique_candidate <- paste0(candidate, "_", suffix)
      suffix <- suffix + 1L
    }
    renamed_id[i] <- unique_candidate
    used_ids <- c(used_ids, unique_candidate)
  }

  source_pairs[, "renamed_id" := renamed_id]
  source_pairs[, "first_source_index" := NULL]
  source_pairs[]
}


apply_feature_merge_rename_map <- function(dt,
                                           rename_map,
                                           id_cols = c(
                                             "feature_id",
                                             "name",
                                             "gene_name",
                                             "exon_id",
                                             "gene_id",
                                             "transcript_id",
                                             "parent_id"
                                           )) {
  if (is.null(dt) || nrow(dt) == 0L || nrow(rename_map) == 0L) {
    return(dt)
  }

  lookup_key <- paste(
    rename_map[["merge_source_index"]],
    rename_map[["merge_id"]],
    sep = "\r"
  )
  lookup <- stats::setNames(rename_map[["renamed_id"]], lookup_key)
  changed_row <- rep(FALSE, nrow(dt))

  for (id_col in intersect(id_cols, names(dt))) {
    value <- as.character(dt[[id_col]])
    key <- paste(dt[[".merge_source_index"]], value, sep = "\r")
    replacement <- unname(lookup[key])
    replace_index <- !is.na(replacement)
    if (any(replace_index)) {
      data.table::set(
        dt,
        i = which(replace_index),
        j = id_col,
        value = replacement[replace_index]
      )
      changed_row[replace_index] <- TRUE
    }
  }

  # Existing raw GFF/GTF attribute strings would contain stale identifiers.
  # Clearing them lets write_feature() reconstruct attributes from normalized
  # identifier columns.
  if ("attribute" %in% names(dt) && any(changed_row)) {
    data.table::set(
      dt,
      i = which(changed_row),
      j = "attribute",
      value = NA_character_
    )
  }
  if ("exon_id" %in% names(dt) && any(changed_row)) {
    data.table::set(
      dt,
      i = which(changed_row),
      j = "exon_id",
      value = NA_character_
    )
  }
  dt
}


sort_feature_merge_table <- function(dt, columns) {
  if (is.null(dt) || nrow(dt) == 0L) {
    return(dt)
  }
  order_columns <- intersect(columns, names(dt))
  if (length(order_columns) > 0L) {
    data.table::setorderv(dt, order_columns)
  }
  dt
}


strip_feature_merge_columns <- function(dt) {
  if (!is.null(dt) && ".merge_source_index" %in% names(dt)) {
    dt[, ".merge_source_index" := NULL]
  }
  dt
}


#' Merge annotation feature objects
#'
#' @description
#' `merge_feature()` combines annotations read from GenePred, GTF, GFF, or BED
#' files into one unified `Feature` object for downstream retrieval and plotting.
#' Inputs may represent different genes, genomic regions, annotation formats, or
#' independently retrieved subsets.
#'
#' Duplicate identifiers are detected across input objects. For non-error
#' conflict strategies, a warning reports duplicated gene, transcript, and
#' feature IDs before the selected strategy is applied. With
#' `conflict = "error"`, the function stops directly without a preliminary
#' warning.
#'
#' With `conflict = "deduplicate"`, input order defines precedence. The first
#' complete copy of a duplicated transcript is retained, duplicated feature
#' records are removed, and gene ranges are recalculated from all retained
#' non-duplicated transcripts. With `conflict = "rename"`, only conflicting
#' identifiers in later inputs are renamed, and all hierarchical references are
#' updated together.
#'
#' @param ... One or more `Feature`, `FeatureTrack`, or `GenePred` objects.
#' A single list containing such objects is also accepted.
#' @param source_names Optional unique source labels. If omitted, names from a
#' supplied list or named arguments are used where available; remaining labels
#' are generated as `track1`, `track2`, and so on. The labels are stored in the
#' `track_source` column.
#' @param conflict Duplicate-ID strategy. `"deduplicate"` keeps the first
#' compatible annotation model according to input order. `"rename"` renames
#' conflicting IDs in later inputs. `"keep_all"` retains conflicts unchanged,
#' and `"error"` stops after duplicate detection. `"keep_first"` is accepted as
#' a backward-compatible alias of `"deduplicate"`.
#' @param rename_prefix Optional character vector with one prefix per input.
#' Used only when `conflict = "rename"`. By default, `source_names` followed by
#' an underscore are used.
#' @param sort Logical. Whether to sort the merged feature and hierarchy tables.
#'
#' @return A merged `Feature` object. If transcript and exon tables are present,
#' the result also inherits from `GenePred`.
#'
#' @examples
#' annotation <- read_genepred(
#'   system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR"),
#'   format = "genePredExt", verbose = FALSE, progress = FALSE
#' )
#' gene_a <- retrieve_feature(annotation, gene_id = "GeneA")
#' gene_b <- retrieve_feature(annotation, gene_id = "GeneB")
#' merged_genes <- merge_feature(gene_a, gene_b)
#' merged_genes
#'
#' @export
merge_feature <- function(...,
                          source_names = NULL,
                          conflict = c(
                            "deduplicate",
                            "rename",
                            "keep_all",
                            "error",
                            "keep_first"
                          ),
                          rename_prefix = NULL,
                          sort = TRUE) {
  objects <- normalize_feature_merge_inputs(list(...))
  stop_if_not(
    length(objects) > 0L,
    "At least one Feature-compatible object is required."
  )
  stop_if_not(
    all(vapply(objects, is_feature_merge_object, logical(1L))),
    paste0(
      "All inputs must be Feature, FeatureTrack, or GenePred objects. ",
      "Read annotation files with read_genepred(), read_gtf(), read_gff(), ",
      "or read_bed() before merging."
    )
  )

  conflict <- match.arg(conflict)
  if (identical(conflict, "keep_first")) {
    conflict <- "deduplicate"
  }

  input_names <- names(objects)
  if (is.null(source_names)) {
    source_names <- paste0("track", seq_along(objects))
    if (!is.null(input_names)) {
      named_input <- !is.na(input_names) & nzchar(input_names)
      source_names[named_input] <- input_names[named_input]
    }
  }
  source_names <- as.character(source_names)
  stop_if_not(
    length(source_names) == length(objects),
    "`source_names` must match the number of input objects."
  )
  stop_if_not(
    all(!is.na(source_names) & nzchar(source_names)),
    "`source_names` must contain non-empty values."
  )
  stop_if_not(
    !anyDuplicated(source_names),
    "`source_names` must be unique."
  )

  if (!is.null(rename_prefix)) {
    stop_if_not(
      identical(conflict, "rename"),
      "`rename_prefix` can only be used with `conflict = \"rename\"`."
    )
    rename_prefix <- as.character(rename_prefix)
    stop_if_not(
      length(rename_prefix) == length(objects),
      "`rename_prefix` must match the number of input objects."
    )
    stop_if_not(
      all(!is.na(rename_prefix) & nzchar(rename_prefix)),
      "`rename_prefix` must contain non-empty values."
    )
  } else if (identical(conflict, "rename")) {
    rename_prefix <- paste0(source_names, "_")
  }

  data_list <- vector("list", length(objects))
  gene_list <- vector("list", length(objects))
  transcript_list <- vector("list", length(objects))
  exon_list <- vector("list", length(objects))

  for (i in seq_along(objects)) {
    object <- objects[[i]]

    data <- data.table::copy(as_feature_table(object))
    data[, "track_source" := source_names[i]]
    data[, ".merge_source_index" := as.integer(i)]
    data_list[[i]] <- data

    genes <- data.table::copy(data.table::as.data.table(object$genes))
    if (nrow(genes) > 0L) {
      genes[, "track_source" := source_names[i]]
      genes[, ".merge_source_index" := as.integer(i)]
    }
    gene_list[[i]] <- genes

    transcripts <- data.table::copy(
      data.table::as.data.table(object$transcripts)
    )
    if (nrow(transcripts) > 0L) {
      transcripts[, "track_source" := source_names[i]]
      transcripts[, ".merge_source_index" := as.integer(i)]
    }
    transcript_list[[i]] <- transcripts

    exons <- data.table::copy(data.table::as.data.table(object$exons))
    if (nrow(exons) > 0L) {
      exons[, "track_source" := source_names[i]]
      exons[, ".merge_source_index" := as.integer(i)]
    }
    exon_list[[i]] <- exons
  }

  data <- data.table::rbindlist(data_list, use.names = TRUE, fill = TRUE)
  genes <- data.table::rbindlist(gene_list, use.names = TRUE, fill = TRUE)
  transcripts <- data.table::rbindlist(
    transcript_list,
    use.names = TRUE,
    fill = TRUE
  )
  exons <- data.table::rbindlist(exon_list, use.names = TRUE, fill = TRUE)

  gene_ownership <- collect_merge_id_sources(
    list(genes, data),
    "gene_id"
  )
  transcript_ownership <- collect_merge_id_sources(
    list(transcripts, exons, data),
    "transcript_id"
  )
  feature_ownership <- collect_merge_id_sources(
    list(data),
    "feature_id"
  )

  all_id_sources <- collect_merge_id_sources(
    list(
      make_merge_id_source_table(genes, "gene_id"),
      make_merge_id_source_table(transcripts, "gene_id"),
      make_merge_id_source_table(exons, "gene_id"),
      make_merge_id_source_table(data, "gene_id"),
      make_merge_id_source_table(transcripts, "transcript_id"),
      make_merge_id_source_table(exons, "transcript_id"),
      make_merge_id_source_table(data, "transcript_id"),
      make_merge_id_source_table(data, "feature_id")
    ),
    "merge_id"
  )

  source_pairs <- data.table::rbindlist(
    list(
      make_merge_id_source_table(genes, "gene_id"),
      make_merge_id_source_table(transcripts, "gene_id"),
      make_merge_id_source_table(exons, "gene_id"),
      make_merge_id_source_table(data, "gene_id"),
      make_merge_id_source_table(transcripts, "transcript_id"),
      make_merge_id_source_table(exons, "transcript_id"),
      make_merge_id_source_table(data, "transcript_id"),
      make_merge_id_source_table(data, "feature_id")
    ),
    use.names = TRUE,
    fill = TRUE
  )
  data.table::setnames(
    source_pairs,
    ".merge_source_index",
    "merge_source_index"
  )
  source_pairs <- unique(
    source_pairs[nonempty_merge_id(merge_id)]
  )
  attr(all_id_sources, "source_pairs") <- source_pairs

  duplicate_ids <- list(
    gene_id = gene_ownership[source_count > 1L, merge_id],
    transcript_id = transcript_ownership[source_count > 1L, merge_id],
    feature_id = feature_ownership[source_count > 1L, merge_id]
  )
  has_duplicates <- any(lengths(duplicate_ids) > 0L)

  if (has_duplicates && identical(conflict, "error")) {
    stop(
      "Duplicated annotation identifiers cannot be merged under strict conflict handling.",
      call. = FALSE
    )
  }

  if (has_duplicates) {
    warning(
      format_merge_duplicate_warning(duplicate_ids, conflict),
      call. = FALSE
    )
  }

  if (has_duplicates && identical(conflict, "deduplicate")) {
    deduplicated <- deduplicate_feature_merge(
      data = data,
      genes = genes,
      transcripts = transcripts,
      exons = exons,
      gene_ownership = gene_ownership,
      transcript_ownership = transcript_ownership,
      feature_ownership = feature_ownership
    )
    data <- deduplicated$data
    genes <- deduplicated$genes
    transcripts <- deduplicated$transcripts
    exons <- deduplicated$exons
  }

  if (has_duplicates && identical(conflict, "rename")) {
    all_existing_ids <- unique(c(
      as.character(genes[["gene_id"]]),
      as.character(transcripts[["transcript_id"]]),
      as.character(data[["feature_id"]])
    ))
    rename_map <- make_feature_merge_rename_map(
      id_sources = all_id_sources,
      rename_prefix = rename_prefix,
      all_existing_ids = all_existing_ids
    )

    data <- apply_feature_merge_rename_map(data, rename_map)
    genes <- apply_feature_merge_rename_map(genes, rename_map)
    transcripts <- apply_feature_merge_rename_map(transcripts, rename_map)
    exons <- apply_feature_merge_rename_map(exons, rename_map)
  }

  if (isTRUE(sort)) {
    data <- sort_feature_merge_table(
      data,
      c(
        "chrom",
        "start",
        "end",
        "gene_id",
        "transcript_id",
        "feature_id",
        "track_source"
      )
    )
    genes <- sort_feature_merge_table(
      genes,
      c("chrom", "gene_start", "gene_end", "gene_id")
    )
    transcripts <- sort_feature_merge_table(
      transcripts,
      c(
        "chrom",
        "tx_start",
        "tx_end",
        "gene_id",
        "transcript_id"
      )
    )
    exons <- sort_feature_merge_table(
      exons,
      c(
        "chrom",
        "exon_start",
        "exon_end",
        "gene_id",
        "transcript_id",
        "exon_number"
      )
    )
  }

  data <- strip_feature_merge_columns(data)
  genes <- strip_feature_merge_columns(genes)
  transcripts <- strip_feature_merge_columns(transcripts)
  exons <- strip_feature_merge_columns(exons)

  source_meta <- lapply(objects, function(object) object$meta %||% list())
  names(source_meta) <- source_names

  output <- Feature(
    data = data,
    genes = genes,
    transcripts = transcripts,
    exons = exons,
    meta = list(
      format = "merged",
      coordinate_internal = "1-based closed",
      source_names = source_names,
      source_meta = source_meta,
      conflict = conflict,
      duplicated_ids = duplicate_ids
    ),
    validation = make_empty_validation()
  )

  if (nrow(transcripts) > 0L && nrow(exons) > 0L) {
    class(output) <- unique(c(class(output), "GenePred"))
  }
  output
}
