# Author: Rensc
# Date: 2026-11-21
# Version: dev001
# Function: Project flat exon and CDS records into transcript groups.
# Input: Feature-compatible annotations with explicit member identifiers.
# Output: A GRangesList with source-row indices for each group.

transcript_view_ids <- function(values) {
  values <- as.character(values)
  values[values %in% c("", ".")] <- NA_character_
  values
}

gff_transcript_key <- function(ids) gsub("(%[[:xdigit:]]{2})", "\\U\\1", ids, perl = TRUE)

feature_transcript_parents <- function(metadata, member_rows, parent_targets,
                                       transcript_ids, parent_key, from_gff) {
  parents <- metadata[["parent_ids"]]
  if (is.null(parents)) {
    parents <- vector("list", length(member_rows))
  } else {
    if (inherits(parents, "List")) parents <- as.list(parents)
    stop_if_not(is.list(parents), "`parent_ids` must be a character-list column.")
    parents <- lapply(parents[member_rows], function(values) {
      if (is.null(values)) return(character())
      stop_if_not(is.character(values) && !anyNA(transcript_view_ids(values)),
        "`parent_ids` must contain nonempty identifier tokens.")
      unique(parent_key(values))
    })
  }

  # Resolve GFF3 Parent tokens against feature IDs before applying display IDs.
  explicit <- lengths(parents) > 0L
  tokens <- unlist(parents[explicit], use.names = FALSE)
  parent_index <- match(tokens, names(parent_targets))
  unknown <- unique(tokens[is.na(parent_index)])
  stop_if_not(!length(unknown),
    sprintf("Unknown transcript Parent ID: %s.", paste(unknown, collapse = ", ")))
  parents[explicit] <- utils::relist(unname(parent_targets[parent_index]), parents[explicit])

  # Legacy annotations and GTF files can explicitly identify transcript members
  # even when their transcript declaration rows are absent.
  legacy_rows <- member_rows[!explicit]
  raw_parent <- transcript_view_ids(metadata$parent_id[legacy_rows])
  resolved_parent <- unname(parent_targets[parent_key(raw_parent)])
  legacy_transcript <- transcript_view_ids(metadata$transcript_id[legacy_rows])
  single_parent <- if (from_gff) {
    pick_first_nonempty(resolved_parent, raw_parent, legacy_transcript)
  } else {
    pick_first_nonempty(legacy_transcript, resolved_parent, raw_parent)
  }
  resolved_single <- unname(parent_targets[parent_key(single_parent)])
  use_raw <- !is.na(resolved_single) & !single_parent %in% transcript_ids
  single_parent[use_raw] <- resolved_single[use_raw]
  missing <- is.na(single_parent)
  stop_if_not(!any(missing), sprintf("Missing transcript assignment for exon or CDS: %s.",
    paste(metadata$feature_id[legacy_rows[missing]], collapse = ", ")))
  parents[!explicit] <- as.list(single_parent)
  parents
}

validate_transcript_group_locations <- function(ranges, rows, ids) {
  locations <- data.table::data.table(transcript = ids,
    chrom = as.character(GenomicRanges::seqnames(ranges))[rows],
    strand = as.character(GenomicRanges::strand(ranges))[rows])
  checked <- locations[, .(conflict = data.table::uniqueN(chrom) > 1L ||
    data.table::uniqueN(strand[strand != "*"]) > 1L), by = "transcript"]
  stop_if_not(!any(checked$conflict), sprintf(
    "Transcript assignments use an incompatible chromosome or strand: %s.",
    paste(checked$transcript[checked$conflict], collapse = ", ")))
}

feature_as_transcript_granges <- function(object) {
  if (inherits(object, "GenePred") && is.null(object$data)) object <- as_feature.GenePred(object)
  ranges <- feature_as_granges(object)
  metadata <- S4Vectors::mcols(ranges)
  format <- toupper(as.character(object$meta$format)[1L])
  from_gff <- format %in% c("GFF", "GFF3")
  parent_key <- if (from_gff || !is.null(metadata$parent_ids)) gff_transcript_key else identity
  types <- tolower(as.character(metadata$type %||% character(length(ranges))))
  feature_levels <- tolower(as.character(metadata$level %||% infer_feature_level(types)))
  transcript_types <- c("snrna", "snorna", "scrna", "tmrna", "lncrna",
    "pseudogenic_transcript", "so:0000234", "so:0000673")
  tx_rows <- which(feature_levels == "transcript" | is_transcript_feature_type(types) |
    types %in% transcript_types)
  member_rows <- which(types %in% c("exon", "cds", "so:0000147", "so:0000316") |
    feature_levels == "exon")

  declared_ids <- transcript_view_ids(metadata$feature_id[tx_rows])
  tx_ids <- pick_first_nonempty(transcript_view_ids(metadata$transcript_id[tx_rows]), declared_ids)
  stop_if_not(!anyNA(tx_ids), "A transcript declaration is missing an identifier.")
  # A discontinuous transcript can have several rows with the same identifiers.
  declarations <- unique(data.frame(raw_id = parent_key(declared_ids), transcript = tx_ids))
  raw <- declarations$raw_id
  duplicate_ids <- unique(c(declarations$transcript[duplicated(declarations$transcript)],
    raw[!is.na(raw) & duplicated(raw)]))
  stop_if_not(!length(duplicate_ids), sprintf("Ambiguous transcript ID: %s.",
    paste(duplicate_ids, collapse = ", ")))
  declared <- !is.na(raw)
  parent_targets <- stats::setNames(declarations$transcript[declared], raw[declared])
  parents <- feature_transcript_parents(metadata, member_rows, parent_targets,
    tx_ids, parent_key, from_gff)
  source_rows <- rep.int(member_rows, lengths(parents))
  members <- as.character(unlist(parents, use.names = FALSE))
  validate_transcript_group_locations(ranges, c(tx_rows, source_rows), c(tx_ids, members))

  groups <- factor(members, levels = unique(c(tx_ids, members)))
  view <- GenomicRanges::split(ranges[source_rows], groups, drop = FALSE)
  S4Vectors::mcols(view) <- S4Vectors::DataFrame(
    source_rows = IRanges::IntegerList(split(source_rows, groups, drop = FALSE))
  )
  S4Vectors::metadata(view) <- S4Vectors::metadata(ranges)
  view
}
