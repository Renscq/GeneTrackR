# Author: Rensc
# Date: 2026-08-31
# Version: dev004
# Function: Extract and plot unified gene feature length distributions
# Input: GenePred object and feature selection
# Output: Length table or ggplot object

#' Extract length records from a gene annotation object
#'
#' @description
#' Build a tidy feature-length table from a Feature/GenePred-compatible annotation object. The function supports
#' gene spans, transcript lengths, exon lengths, CDS lengths, total UTR lengths,
#' 5' UTR lengths, and 3' UTR lengths. CDS and UTR lengths can be summarized per
#' transcript or returned as individual genomic segments.
#'
#' @param object A Feature or GenePred-compatible annotation object.
#' @param feature Feature type to extract. Use `all` to return gene, transcript,
#' exon, CDS, UTR, 5' UTR, and 3' UTR records.
#' @param unit Output unit. `auto` uses gene-level records for genes,
#' transcript-level records for transcripts, transcript-level totals for CDS/UTR,
#' and segment-level records for exons. `transcript` summarizes exonic features
#' per transcript where applicable. `segment` returns individual exon/CDS/UTR
#' genomic segments where applicable.
#' @param transcript_length Transcript length definition. `spliced` uses the sum
#' of exon lengths. `genomic` uses transcript span length.
#' @param chrom Optional chromosome filter.
#' @param start Optional region start in 1-based closed coordinates.
#' @param end Optional region end in 1-based closed coordinates.
#' @param mode Region selection mode passed to [retrieve_feature()].
#' @param keep_zero Logical. Whether to keep zero-length transcript-level CDS/UTR
#' records. Default is `FALSE`.
#' @param ... Deprecated compatibility arguments forwarded by
#' `get_genepred_length_table()` to `get_gene_length_distribution_table()`.
#' @details
#' Use `unit = "segment"` to inspect each exon/CDS/UTR segment separately,
#' and `unit = "transcript"` to summarize feature lengths per transcript.
#' `transcript_length = "spliced"` uses exon lengths, whereas `"genomic"`
#' uses transcript span length including introns.
#' @return A data.table with feature-length records.
#' @examples
#' gp <- read_genepred(
#'   system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR"),
#'   format = "genePredExt", verbose = FALSE, progress = FALSE
#' )
#' get_gene_length_distribution_table(gp, feature = "cds", unit = "transcript")
#' @export
get_gene_length_distribution_table <- function(object,
                                      feature = c("all", "gene", "transcript", "exon", "cds", "utr", "five_utr", "three_utr"),
                                      unit = c("auto", "transcript", "segment"),
                                      transcript_length = c("spliced", "genomic"),
                                      chrom = NULL,
                                      start = NULL,
                                      end = NULL,
                                      mode = c("overlap", "within", "trim"),
                                      keep_zero = FALSE) {
  stop_if_not(is_gene_model_feature(object), "`object` must be a Feature or GenePred-compatible annotation object with transcript/exon records.")
  object <- as_genepred(object)

  feature <- match.arg(feature)
  unit <- match.arg(unit)
  transcript_length <- match.arg(transcript_length)
  mode <- match.arg(mode)
  check_region(chrom, start, end)

  obj <- if (!is.null(chrom)) {
    retrieve_feature(
      object = object,
      chrom = chrom,
      start = start %||% 1L,
      end = end %||% .Machine$integer.max,
      mode = mode,
      as = "Feature"
    )
  } else {
    object
  }

  selected_features <- if (feature == "all") {
    c("gene", "transcript", "exon", "cds", "utr", "five_utr", "three_utr")
  } else {
    feature
  }

  out <- lapply(selected_features, function(x) {
    get_single_gene_length_distribution_table(
      object = obj,
      feature = x,
      unit = unit,
      transcript_length = transcript_length,
      keep_zero = keep_zero
    )
  })

  out <- data.table::rbindlist(out, fill = TRUE)
  if (nrow(out) == 0L) {
    return(empty_gene_length_distribution_table())
  }

  out[, feature := factor(
    feature,
    levels = c("gene", "transcript", "exon", "cds", "utr", "five_utr", "three_utr")
  )]
  data.table::setorder(out, feature, chrom, start, end, gene_id, transcript_id)
  out[]
}

#' Plot gene annotation feature length distributions
#'
#' @description
#' Plot the length distribution of genes, transcripts, exons, CDS, UTRs, 5' UTRs,
#' and 3' UTRs from any Feature/GenePred-compatible annotation object. The function is designed for annotation
#' quality control and comparison of coding versus non-coding feature lengths.
#'
#' @param object A Feature or GenePred-compatible annotation object.
#' @param feature Feature type to plot. Use `all` to facet multiple feature types.
#' @param unit Output unit used for length extraction. See
#' [get_gene_length_distribution_table()].
#' @param transcript_length Transcript length definition for transcript records.
#' @param chrom Optional chromosome filter.
#' @param start Optional region start in 1-based closed coordinates.
#' @param end Optional region end in 1-based closed coordinates.
#' @param mode Region selection mode passed to [retrieve_feature()].
#' @param group_by Grouping variable. Common choices are `gene_type`, `feature`,
#' `strand`, and `chrom`.
#' @param plot_type Plot type. One of `density`, `histogram`, `boxplot`, or
#' `violin`.
#' @param scale Length scale. `log10` is recommended for genomic features because
#' feature lengths are usually right-skewed.
#' @param bins Number of bins for histogram.
#' @param facet Logical. Whether to facet by feature when multiple features are
#' requested.
#' @param keep_zero Logical. Whether to keep zero-length CDS/UTR records.
#' @param fill_palette RColorBrewer palette name used for grouped fills. Default is `Paired`.
#' @param fill_colors Optional custom fill colors. Named vectors are matched to the
#' values of `group_by`, for example `c(coding = "#1f78b4", `non-coding` = "#a6cee3")`. Unnamed colors are matched in plotting order and are automatically extended when needed.
#' @param border_color Optional border color for histograms, boxplots, and violin plots.
#' Use `NA` to hide borders.
#' @param return_data Logical. If `TRUE`, return a list containing the plot and
#' the underlying length table.
#' @param ... Deprecated compatibility arguments forwarded by
#' `plot_genepred_length_distribution()` to `plot_gene_length_distribution()`.
#' @details
#' `fill_colors` follows the values of `group_by`. For example, if
#' `group_by = "gene_type"`, names should match `coding` and `non-coding`; if
#' `group_by = "feature"`, names can include `gene`, `transcript`, `exon`,
#' `cds`, `utr`, `five_utr`, and `three_utr`. Use `return_data = TRUE` to inspect
#' the exact groups before assigning named colors. See also
#' [GeneTrackR-advanced-parameters].
#' @return A ggplot object, or a list with `plot` and `data` when
#' `return_data = TRUE`.
#' @examples
#' gp <- read_genepred(
#'   system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR"),
#'   format = "genePredExt", verbose = FALSE, progress = FALSE
#' )
#' plot_gene_length_distribution(gp, feature = "gene", plot_type = "histogram")
#' @export
plot_gene_length_distribution <- function(object,
                                              feature = c("all", "gene", "transcript", "exon", "cds", "utr", "five_utr", "three_utr"),
                                              unit = c("auto", "transcript", "segment"),
                                              transcript_length = c("spliced", "genomic"),
                                              chrom = NULL,
                                              start = NULL,
                                              end = NULL,
                                              mode = c("overlap", "within", "trim"),
                                              group_by = c("gene_type", "feature", "strand", "chrom", "none"),
                                              plot_type = c("density", "histogram", "boxplot", "violin"),
                                              scale = c("log10", "linear"),
                                              bins = 60L,
                                              facet = TRUE,
                                              keep_zero = FALSE,
                                              fill_palette = "Paired",
                                              fill_colors = NULL,
                                              border_color = NA,
                                              return_data = FALSE) {
  feature <- match.arg(feature)
  unit <- match.arg(unit)
  transcript_length <- match.arg(transcript_length)
  mode <- match.arg(mode)
  group_by <- match.arg(group_by)
  plot_type <- match.arg(plot_type)
  scale <- match.arg(scale)

  dt <- get_gene_length_distribution_table(
    object = object,
    feature = feature,
    unit = unit,
    transcript_length = transcript_length,
    chrom = chrom,
    start = start,
    end = end,
    mode = mode,
    keep_zero = keep_zero
  )

  if (nrow(dt) == 0L) {
    stop("No length records are available for the requested feature and region.", call. = FALSE)
  }

  dt <- dt[length_bp > 0]
  if (nrow(dt) == 0L) {
    stop("No positive-length records are available for plotting.", call. = FALSE)
  }

  dt[, length_plot := if (scale == "log10") log10(length_bp) else as.numeric(length_bp)]
  x_label <- if (scale == "log10") "log10(length, bp)" else "Length (bp)"

  if (group_by == "none") {
    dt[, plot_group := "all"]
    group_col <- "plot_group"
  } else {
    group_col <- group_by
  }

  p <- build_gene_length_distribution_plot(
    dt = dt,
    group_col = group_col,
    plot_type = plot_type,
    bins = bins,
    x_label = x_label,
    fill_palette = fill_palette,
    fill_colors = fill_colors,
    border_color = border_color
  )

  if (facet && length(unique(as.character(dt$feature))) > 1L && plot_type %in% c("density", "histogram")) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$feature), scales = "free_y")
  }

  p <- p +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      legend.title = ggplot2::element_blank()
    )

  if (return_data) {
    return(list(plot = p, data = dt))
  }
  p
}


#' @rdname get_gene_length_distribution_table
#' @export
get_genepred_length_table <- function(...) {
  get_gene_length_distribution_table(...)
}

#' @rdname plot_gene_length_distribution
#' @export
plot_genepred_length_distribution <- function(...) {
  plot_gene_length_distribution(...)
}

get_single_gene_length_distribution_table <- function(object, feature, unit, transcript_length, keep_zero) {
  if (feature == "gene") {
    return(get_gene_length_table(object))
  }

  if (feature == "transcript") {
    return(get_transcript_length_table(object, transcript_length = transcript_length))
  }

  if (feature == "exon") {
    if (unit == "transcript") {
      return(get_exon_total_length_table(object))
    }
    return(get_exon_segment_length_table(object))
  }

  if (feature %in% c("cds", "utr", "five_utr", "three_utr")) {
    segment_dt <- get_exonic_subfeature_segments(object, feature = feature)
    if (unit == "segment") {
      return(filter_zero_length(segment_dt, keep_zero = keep_zero))
    }
    return(summarise_segments_by_transcript(segment_dt, feature = feature, keep_zero = keep_zero))
  }

  empty_gene_length_distribution_table()
}

get_gene_length_table <- function(object) {
  dt <- as_gene_table(object)
  if (nrow(dt) == 0L) {
    return(empty_gene_length_distribution_table())
  }

  dt[, .(
    feature = "gene",
    unit = "gene",
    feature_id = gene_id,
    gene_id = gene_id,
    transcript_id = NA_character_,
    chrom = chrom,
    strand = strand,
    gene_type = gene_type,
    start = as.integer(gene_start),
    end = as.integer(gene_end),
    length_bp = as.numeric(gene_end - gene_start + 1L)
  )]
}

get_transcript_length_table <- function(object, transcript_length = c("spliced", "genomic")) {
  transcript_length <- match.arg(transcript_length)
  tx <- as_transcript_table(object)
  ex <- as_exon_table(object)
  if (nrow(tx) == 0L) {
    return(empty_gene_length_distribution_table())
  }

  if (transcript_length == "spliced") {
    len <- ex[, .(
      length_bp = as.numeric(sum(as.numeric(exon_end - exon_start + 1L), na.rm = TRUE))
    ), by = transcript_id]
    tx <- merge(tx, len, by = "transcript_id", all.x = TRUE)
  } else {
    tx[, length_bp := as.numeric(tx_end - tx_start + 1L)]
  }

  tx[is.na(length_bp), length_bp := 0]
  tx[, .(
    feature = "transcript",
    unit = "transcript",
    feature_id = transcript_id,
    gene_id = gene_id,
    transcript_id = transcript_id,
    chrom = chrom,
    strand = strand,
    gene_type = gene_type,
    start = as.integer(tx_start),
    end = as.integer(tx_end),
    length_bp = as.numeric(length_bp)
  )]
}

get_exon_segment_length_table <- function(object) {
  ex <- as_exon_table(object)
  tx <- as_transcript_table(object)[, .(transcript_id, gene_type)]
  if (nrow(ex) == 0L) {
    return(empty_gene_length_distribution_table())
  }

  ex <- merge(ex, tx, by = "transcript_id", all.x = TRUE)
  ex[, .(
    feature = "exon",
    unit = "segment",
    feature_id = paste(transcript_id, exon_number, sep = ":exon"),
    gene_id = gene_id,
    transcript_id = transcript_id,
    chrom = chrom,
    strand = strand,
    gene_type = gene_type,
    start = as.integer(exon_start),
    end = as.integer(exon_end),
    length_bp = as.numeric(exon_end - exon_start + 1L)
  )]
}

get_exon_total_length_table <- function(object) {
  ex <- get_exon_segment_length_table(object)
  if (nrow(ex) == 0L) {
    return(empty_gene_length_distribution_table())
  }
  summarise_segments_by_transcript(ex, feature = "exon", keep_zero = FALSE)
}

get_exonic_subfeature_segments <- function(object, feature = c("cds", "utr", "five_utr", "three_utr")) {
  feature <- match.arg(feature)
  tx <- as_transcript_table(object)
  ex <- as_exon_table(object)

  if (nrow(tx) == 0L || nrow(ex) == 0L) {
    return(empty_gene_length_distribution_table())
  }

  tx <- tx[, .(transcript_id, cds_start, cds_end, gene_type)]
  ex <- merge(ex, tx, by = "transcript_id", all.x = TRUE)
  ex <- ex[gene_type == "coding" & !is.na(cds_start) & !is.na(cds_end) & cds_start <= cds_end]

  if (nrow(ex) == 0L) {
    return(empty_gene_length_distribution_table())
  }

  if (feature == "cds") {
    seg <- make_interval_overlap_segments(
      ex,
      feature = "cds",
      interval_start = ex$cds_start,
      interval_end = ex$cds_end
    )
    return(seg)
  }

  five <- data.table::rbindlist(list(
    make_interval_overlap_segments(
      ex[strand == "+"],
      feature = "five_utr",
      interval_start = -Inf,
      interval_end = ex[strand == "+"]$cds_start - 1L
    ),
    make_interval_overlap_segments(
      ex[strand == "-"],
      feature = "five_utr",
      interval_start = ex[strand == "-"]$cds_end + 1L,
      interval_end = Inf
    )
  ), fill = TRUE)

  three <- data.table::rbindlist(list(
    make_interval_overlap_segments(
      ex[strand == "+"],
      feature = "three_utr",
      interval_start = ex[strand == "+"]$cds_end + 1L,
      interval_end = Inf
    ),
    make_interval_overlap_segments(
      ex[strand == "-"],
      feature = "three_utr",
      interval_start = -Inf,
      interval_end = ex[strand == "-"]$cds_start - 1L
    )
  ), fill = TRUE)

  if (feature == "five_utr") {
    return(five)
  }
  if (feature == "three_utr") {
    return(three)
  }

  utr <- data.table::rbindlist(list(five, three), fill = TRUE)
  if (nrow(utr) > 0L) {
    utr[, feature := "utr"]
    utr[, feature_id := paste(transcript_id, seq_len(.N), sep = ":utr"), by = transcript_id]
  }
  utr
}

make_interval_overlap_segments <- function(ex, feature, interval_start, interval_end) {
  if (nrow(ex) == 0L) {
    return(empty_gene_length_distribution_table())
  }

  seg <- data.table::copy(ex)
  seg[, start := as.integer(pmax(as.numeric(exon_start), as.numeric(interval_start)))]
  seg[, end := as.integer(pmin(as.numeric(exon_end), as.numeric(interval_end)))]
  seg <- seg[start <= end]

  if (nrow(seg) == 0L) {
    return(empty_gene_length_distribution_table())
  }

  seg[, .(
    feature = feature,
    unit = "segment",
    feature_id = paste(transcript_id, exon_number, feature, sep = ":"),
    gene_id = gene_id,
    transcript_id = transcript_id,
    chrom = chrom,
    strand = strand,
    gene_type = gene_type,
    start = as.integer(start),
    end = as.integer(end),
    length_bp = as.numeric(end - start + 1L)
  )]
}

summarise_segments_by_transcript <- function(segment_dt, feature, keep_zero = FALSE) {
  if (nrow(segment_dt) == 0L) {
    return(empty_gene_length_distribution_table())
  }

  out <- segment_dt[, .(
    feature = feature,
    unit = "transcript",
    feature_id = transcript_id[1L],
    gene_id = gene_id[1L],
    transcript_id = transcript_id[1L],
    chrom = chrom[1L],
    strand = strand[1L],
    gene_type = gene_type[1L],
    start = as.integer(min(start, na.rm = TRUE)),
    end = as.integer(max(end, na.rm = TRUE)),
    length_bp = as.numeric(sum(length_bp, na.rm = TRUE))
  ), by = transcript_id]

  out[, transcript_id := NULL]
  data.table::setnames(out, "transcript_id", "group_transcript_id", skip_absent = TRUE)

  out[, transcript_id := feature_id]
  out <- out[, .(
    feature,
    unit,
    feature_id,
    gene_id,
    transcript_id,
    chrom,
    strand,
    gene_type,
    start,
    end,
    length_bp
  )]

  filter_zero_length(out, keep_zero = keep_zero)
}

filter_zero_length <- function(dt, keep_zero = FALSE) {
  if (keep_zero || nrow(dt) == 0L) {
    return(dt)
  }
  dt[length_bp > 0]
}

empty_gene_length_distribution_table <- function() {
  data.table::data.table(
    feature = character(),
    unit = character(),
    feature_id = character(),
    gene_id = character(),
    transcript_id = character(),
    chrom = character(),
    strand = character(),
    gene_type = character(),
    start = integer(),
    end = integer(),
    length_bp = numeric()
  )
}

build_gene_length_distribution_plot <- function(dt, group_col, plot_type, bins, x_label, fill_palette = "Paired", fill_colors = NULL, border_color = NA) {
  border_color <- normalize_border_color(border_color)

  if (plot_type == "density") {
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = .data$length_plot, fill = .data[[group_col]])) +
      ggplot2::geom_density(alpha = 0.35, color = border_color, linewidth = 0.6, na.rm = TRUE) +
      ggplot2::labs(x = x_label, y = "Density")
    return(apply_discrete_fill_scale(p, color_palette = fill_palette, fill_colors = fill_colors))
  }

  if (plot_type == "histogram") {
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = .data$length_plot, fill = .data[[group_col]])) +
      ggplot2::geom_histogram(bins = bins, alpha = 0.7, position = "identity", color = border_color, na.rm = TRUE) +
      ggplot2::labs(x = x_label, y = "Count")
    return(apply_discrete_fill_scale(p, color_palette = fill_palette, fill_colors = fill_colors))
  }

  if (plot_type == "boxplot") {
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = .data[[group_col]], y = .data$length_plot, fill = .data[[group_col]])) +
      ggplot2::geom_boxplot(color = border_color, outlier.alpha = 0.25, na.rm = TRUE) +
      ggplot2::labs(x = group_col, y = x_label)
    return(apply_discrete_fill_scale(p, color_palette = fill_palette, fill_colors = fill_colors))
  }

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = .data[[group_col]], y = .data$length_plot, fill = .data[[group_col]])) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.7, color = border_color, na.rm = TRUE) +
    ggplot2::geom_boxplot(width = 0.12, color = border_color, outlier.alpha = 0.15, na.rm = TRUE) +
    ggplot2::labs(x = group_col, y = x_label)
  apply_discrete_fill_scale(p, color_palette = fill_palette, fill_colors = fill_colors)
}
