# Author: Rensc
# Date: 2026-07-30
# Version: dev007
# Function: Plot haplotype-variant and haplotype-phenotype figures
# Input: HapVariant and phenotype objects
# Output: ggplot or patchwork figures

#' Plot gene variants and haplotype table
#'
#' @description
#' Draws a haplotype-variant figure with a compact gene model track, variant
#' markers, connector lines, and a haplotype genotype table. Genomic gaps are
#' compressed by mapping variant positions to evenly spaced table columns while
#' retaining the relative order of gene model features.
#'
#' @param hap A HapVariant object from `hap_variant()`.
#' @param annotation Optional Feature/GenePred annotation object. When supplied,
#' a compact gene track is drawn above the haplotype table.
#' @param show_gene_model Logical. Whether to draw the gene track when
#' `annotation` is supplied.
#' @param min_hap_samples Minimum sample number required for a haplotype group to be displayed.
#' @param show_reference_row Logical. Whether to add two reference rows showing REF and ALT alleles for each variant.
#' @param variant_label Column used for variant labels. One of `variant_id`, `pos`, or an existing column in `hap$variants`.
#' @param show_gene_pos_axis Logical. Whether to show genomic coordinate labels above the gene track.
#' @param gene_pos_axis_n Approximate number of genomic coordinate ticks above the gene track.
#' @param gene_pos_axis_label Optional x-axis title for genomic coordinate labels.
#' @param gene_pos_x_angle Angle of gene-position x-axis labels. Default is 0.
#' @param gene_track_legend_position Legend position for variant-type markers in the gene track. One of `right`, `top`, or `none`.
#' @param direction_mode Gene-strand arrow style for the compact gene track. `transcript` draws one arrow per transcript, `gene` draws one arrow per gene, `end` draws a short arrow at the directional end of each gene, and `none` hides direction arrows.
#' @param text_size Base text size for gene-track labels, axes, legends, and table axes.
#' @param table_x_angle Angle of haplotype table x-axis labels. Default is 90.
#' @param genotype_text_size Genotype cell text size only. Default is 3.2.
#' @param gene_track_height Relative height of the gene track panel.
#' @param connector_height Relative height of the connector panel.
#' @param table_height Relative height of the haplotype table panel.
#' @param exon_height Height of exon boxes in the compact gene track.
#' @param cds_height Height of CDS boxes in the compact gene track.
#' @param gene_palette RColorBrewer palette name used for gene model feature fills.
#' @param gene_colors Optional custom fill colors for gene model features.
#' @param gene_border_color Border color for gene model feature boxes.
#' @param table_palette RColorBrewer palette name used for haplotype table fills.
#' For allele-string haplotypes, the first five colors are always assigned in
#' the fixed order A, T, C, G, and indel.
#' @param table_colors Optional custom fill colors for haplotype table genotypes.
#' For allele-string haplotypes, use a vector named `A`, `T`, `C`, `G`, and
#' `indel`, or supply five unnamed colors in that order.
#' @param table_alpha Alpha value for haplotype table fill colors.
#' @param reference_fill Background fill color for the REF and ALT reference rows.
#' @param variant_palette RColorBrewer palette name used for solid variant-type triangle marker colors.
#' @param variant_colors Optional custom colors for solid variant-type triangle markers.
#' @param variant_alpha Alpha value for solid variant-type triangle marker colors.
#' @param show_variant_marker Logical. Whether to draw natural-variant triangle markers.
#' @param variant_marker_size Size of natural-variant triangle markers. Set to 0 to hide markers.
#' @return A patchwork object with attributes `plot_data`, `variant_data`, and `gene_data`.
#' @examples
#' vcf_file <- system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR")
#' anno_file <- system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file, mode = "memory", verbose = FALSE)
#' anno <- read_genepred(anno_file, format = "genePredExt", verbose = FALSE)
#' hap <- hap_variant(vcf, annotation = anno, gene_id = "GeneA", genotype_mode = "string", min_variant_number = 1)
#' plot_hap_variant(hap, annotation = anno, min_hap_samples = 1)
#' plot_hap_variant(
#'   hap,
#'   annotation = anno,
#'   min_hap_samples = 1,
#'   table_x_angle = 90,
#'   table_palette = "Paired"
#' )
#' @export
plot_hap_variant <- function(hap,
                             annotation = NULL,
                             show_gene_model = TRUE,
                             min_hap_samples = 5L,
                             show_reference_row = TRUE,
                             variant_label = c("variant_id", "pos"),
                             show_gene_pos_axis = TRUE,
                             gene_pos_axis_n = 5L,
                             gene_pos_axis_label = NULL,
                             gene_pos_x_angle = 0,
                             gene_track_legend_position = c("right", "top", "none"),
                             direction_mode = c("transcript", "gene", "end", "none"),
                             text_size = 14,
                             table_x_angle = 90,
                             genotype_text_size = 3.2,
                             gene_track_height = 1.25,
                             connector_height = 0.35,
                             table_height = NULL,
                             exon_height = 0.22,
                             cds_height = 0.44,
                             gene_palette = "Paired",
                             gene_colors = NULL,
                             gene_border_color = NA,
                             table_palette = "Paired",
                             table_colors = NULL,
                             table_alpha = 0.6,
                             reference_fill = "white",
                             variant_palette = "Paired",
                             variant_colors = NULL,
                             variant_alpha = 0.6,
                             show_variant_marker = TRUE,
                             variant_marker_size = 2.8) {
  stop_if_not(inherits(hap, "HapVariant"), "`hap` must be a HapVariant object.")

  vars <- data.table::as.data.table(hap$variants)
  haps <- data.table::as.data.table(hap$haplotypes)
  stop_if_not(nrow(vars) > 0L, "No variants are available in `hap$variants`.")
  stop_if_not(nrow(haps) > 0L, "No haplotypes are available in `hap$haplotypes`.")

  min_hap_samples <- as.integer(min_hap_samples)[1L]
  stop_if_not(!is.na(min_hap_samples) && min_hap_samples >= 1L, "`min_hap_samples` must be a positive integer.")
  show_reference_row <- isTRUE(show_reference_row)
  gene_border_color <- normalize_border_color(gene_border_color)
  table_alpha <- as.numeric(table_alpha)[1L]
  if (is.na(table_alpha)) table_alpha <- 0.6
  table_alpha <- max(0, min(1, table_alpha))
  genotype_text_size <- as.numeric(genotype_text_size)[1L]
  if (is.na(genotype_text_size) || genotype_text_size <= 0) {
    genotype_text_size <- 3.2
  }
  table_x_angle <- as.numeric(table_x_angle)[1L]
  if (is.na(table_x_angle)) table_x_angle <- 90
  table_x_angle <- max(0, min(180, table_x_angle))
  gene_pos_x_angle <- as.numeric(gene_pos_x_angle)[1L]
  if (is.na(gene_pos_x_angle)) gene_pos_x_angle <- 0
  gene_pos_x_angle <- max(0, min(180, gene_pos_x_angle))
  reference_fill <- as.character(reference_fill)[1L]
  if (is.na(reference_fill) || reference_fill == "") {
    reference_fill <- "white"
  }
  variant_alpha <- as.numeric(variant_alpha)[1L]
  if (is.na(variant_alpha)) variant_alpha <- 0.6
  variant_alpha <- max(0, min(1, variant_alpha))
  show_variant_marker <- isTRUE(show_variant_marker)
  variant_marker_size <- as.numeric(variant_marker_size)[1L]
  if (is.na(variant_marker_size)) variant_marker_size <- 2.8
  if (variant_marker_size <= 0) {
    show_variant_marker <- FALSE
    variant_marker_size <- 0
  }
  if (!"variant_type" %in% names(vars)) {
    vars[, "variant_type_label" := "..."]
  } else {
    vars[, "variant_type_label" := normalize_variant_marker_type(variant_type)]
  }

  haps <- haps[sample_n >= min_hap_samples]
  stop_if_not(nrow(haps) > 0L, "No haplotypes remained after `min_hap_samples` filtering.")

  variant_ids <- vars[["variant_id"]]
  variant_ids <- variant_ids[variant_ids %in% names(haps)]
  stop_if_not(length(variant_ids) > 0L, "No variant genotype columns were found in `hap$haplotypes`.")

  vars <- vars[match(variant_ids, variant_id)]
  vars[, "variant_index" := seq_len(.N)]

  variant_label <- match.arg(variant_label)
  gene_track_legend_position <- match.arg(gene_track_legend_position)
  direction_mode <- match.arg(direction_mode)
  show_gene_pos_axis <- isTRUE(show_gene_pos_axis)
  gene_pos_axis_n <- as.integer(gene_pos_axis_n)[1L]
  if (is.na(gene_pos_axis_n) || gene_pos_axis_n < 2L) {
    gene_pos_axis_n <- 5L
  }
  if (is.null(gene_pos_axis_label)) {
    gene_pos_axis_label <- NA_character_
  }
  gene_pos_axis_label <- as.character(gene_pos_axis_label)[1L]

  if (!variant_label %in% names(vars)) {
    variant_label <- "variant_id"
  }
  vars[, "variant_label" := as.character(.SD[[1L]]), .SDcols = variant_label]
  vars[is.na(variant_label) | variant_label == "", "variant_label" := as.character(variant_id)]

  id_vars <- intersect(c("hap_id", "sample_n", "samples"), names(haps))
  table_long <- data.table::melt(
    haps,
    id.vars = id_vars,
    measure.vars = variant_ids,
    variable.name = "variant_id",
    value.name = "genotype",
    variable.factor = FALSE
  )

  table_long <- merge(
    table_long,
    vars[, .(variant_id, variant_index, variant_label, ref, alt)],
    by = "variant_id",
    all.x = TRUE,
    sort = FALSE
  )

  haps[, "hap_sort_number" := parse_hap_numeric_id(hap_id)]
  data.table::setorderv(haps, c("sample_n", "hap_sort_number", "hap_id"), order = c(-1L, 1L, 1L))
  hap_order <- haps$hap_id
  table_long[, "row_label" := paste0(hap_id, " (n=", sample_n, ")")]
  hap_label_levels <- haps[, paste0(hap_id, " (n=", sample_n, ")")]
  row_levels <- rev(hap_label_levels)

  table_long[, "genotype_label" := format_hap_table_genotype_label(
    genotype = genotype,
    ref = ref,
    alt = alt,
    genotype_mode = hap$meta$genotype_mode %||% NA_character_
  )]
  table_long[is.na(genotype_label), "genotype_label" := "NA"]
  table_long[, "genotype_fill_class" := normalize_hap_table_fill_class(
    genotype_label,
    genotype_mode = hap$meta$genotype_mode %||% NA_character_
  )]

  if (show_reference_row) {
    ref_row_label <- "REF"
    alt_row_label <- "ALT"
    ref_row <- vars[, .(
      variant_id,
      variant_index,
      variant_label,
      row_label = ref_row_label,
      genotype = format_hap_allele(ref),
      genotype_label = format_hap_allele(ref),
      allele_label = format_hap_allele(ref)
    )]
    alt_row <- vars[, .(
      variant_id,
      variant_index,
      variant_label,
      row_label = alt_row_label,
      genotype = vapply(strsplit(as.character(alt), ",", fixed = TRUE), function(x) {
        x <- x[!is.na(x) & nzchar(x)]
        x <- format_hap_allele(x)
        x <- x[!is.na(x) & nzchar(x)]
        if (length(x) == 0L) return("NA")
        paste(x, collapse = ",")
      }, character(1L)),
      genotype_label = vapply(strsplit(as.character(alt), ",", fixed = TRUE), function(x) {
        x <- x[!is.na(x) & nzchar(x)]
        x <- format_hap_allele(x)
        x <- x[!is.na(x) & nzchar(x)]
        if (length(x) == 0L) return("NA")
        paste(x, collapse = ",")
      }, character(1L)),
      allele_label = vapply(strsplit(as.character(alt), ",", fixed = TRUE), function(x) {
        x <- x[!is.na(x) & nzchar(x)]
        x <- format_hap_allele(x)
        x <- x[!is.na(x) & nzchar(x)]
        if (length(x) == 0L) return("NA")
        paste(x, collapse = ",")
      }, character(1L))
    )]
    ref_alt <- data.table::rbindlist(list(alt_row, ref_row), fill = TRUE)
    for (nm in setdiff(names(table_long), names(ref_alt))) {
      ref_alt[, (nm) := NA]
    }
    data.table::setcolorder(ref_alt, names(table_long))
    table_long <- data.table::rbindlist(list(table_long, ref_alt), fill = TRUE)
    row_levels <- c(row_levels, alt_row_label, ref_row_label)
  }

  table_long[, "hap_y" := match(row_label, row_levels)]

  x_limits <- c(0.5, length(variant_ids) + 0.5)
  x_breaks <- vars$variant_index
  x_labels <- vars$variant_label

  table_mapper <- make_hap_table_mapper(vars)
  gene_mapper <- make_hap_gene_mapper(vars, hap$region)
  vars[, "table_x" := table_mapper(pos)]
  vars[, "gene_x" := gene_mapper(pos)]

  gene_data <- NULL
  p_gene <- NULL
  if (isTRUE(show_gene_model) && !is.null(annotation)) {
    gene_data <- prepare_hap_gene_track(annotation, hap, vars, gene_mapper)
    if (!is.null(gene_data) && nrow(gene_data$transcripts) > 0L) {
      p_gene <- draw_hap_gene_track(
        gene_data = gene_data,
        vars = vars,
        x_limits = x_limits,
        x_breaks = x_breaks,
        x_labels = x_labels,
        text_size = text_size,
        gene_text_size = text_size,
        exon_height = exon_height,
        cds_height = cds_height,
        gene_palette = gene_palette,
        gene_colors = gene_colors,
        gene_border_color = gene_border_color,
        variant_palette = variant_palette,
        variant_colors = variant_colors,
        variant_alpha = variant_alpha,
        show_variant_marker = show_variant_marker,
        variant_marker_size = variant_marker_size,
        gene_track_legend_position = gene_track_legend_position,
        direction_mode = direction_mode,
        show_gene_pos_axis = show_gene_pos_axis,
        gene_pos_axis_n = gene_pos_axis_n,
        gene_pos_axis_label = gene_pos_axis_label,
        gene_pos_x_angle = gene_pos_x_angle
      )
    }
  }

  p_connector <- ggplot2::ggplot(vars) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$gene_x, xend = .data$table_x, y = 1, yend = 0),
      linewidth = 0.35,
      color = "grey35",
      lineend = "butt"
    ) +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 12, 0, 8))

  n_rows <- length(row_levels)
  table_genotype <- table_long
  table_reference <- table_long[0]
  if (isTRUE(show_reference_row)) {
    table_reference <- table_long[row_label %in% c("REF", "ALT")]
    table_genotype <- table_long[!row_label %in% c("REF", "ALT")]
  }
  if (!"allele_label" %in% names(table_reference)) {
    table_reference[, "allele_label" := character()]
  }

  p_table <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = table_genotype,
      ggplot2::aes(x = .data$variant_index, y = .data$hap_y, fill = .data$genotype_fill_class),
      color = "grey80",
      linewidth = 0.28,
      width = 0.96,
      height = 0.86
    ) +
    ggplot2::geom_tile(
      data = table_reference,
      ggplot2::aes(x = .data$variant_index, y = .data$hap_y),
      fill = reference_fill,
      color = "grey70",
      linewidth = 0.28,
      width = 0.96,
      height = 0.86
    ) +
    ggplot2::geom_text(
      data = table_genotype,
      ggplot2::aes(x = .data$variant_index, y = .data$hap_y, label = .data$genotype_label),
      size = genotype_text_size,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = table_reference,
      ggplot2::aes(x = .data$variant_index, y = .data$hap_y, label = .data$allele_label),
      size = genotype_text_size * 0.95,
      color = "black",
      vjust = 0.5
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = x_limits,
      expand = c(0, 0),
      position = "bottom"
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_len(n_rows),
      labels = row_levels,
      limits = c(0.5, n_rows + 0.5),
      expand = c(0, 0)
    ) +
    ggplot2::labs(x = NULL, y = "Haplotype") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(size = text_size, angle = table_x_angle, hjust = 1, vjust = 0.5, color = "black"),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = text_size, color = "black"),
      axis.title.y = ggplot2::element_text(size = text_size, color = "black"),
      axis.title.x = ggplot2::element_blank(),
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 12, 8, 8)
    )

  p_table <- apply_hap_table_fill_scale(
    p_table,
    values = table_genotype$genotype_label,
    table_palette = table_palette,
    table_colors = table_colors,
    table_alpha = table_alpha,
    fixed_allele_classes = identical(
      as.character(hap$meta$genotype_mode %||% NA_character_)[1L],
      "string"
    )
  )

  if (is.null(table_height)) {
    table_height <- max(1.2, n_rows * 0.38)
  }

  if (!is.null(p_gene)) {
    p <- patchwork::wrap_plots(
      p_gene,
      p_connector,
      p_table,
      ncol = 1,
      heights = c(gene_track_height, connector_height, table_height)
    )
  } else {
    p_variant <- ggplot2::ggplot(vars, ggplot2::aes(x = .data$gene_x, y = 1))
    if (isTRUE(show_variant_marker)) {
      p_variant <- p_variant +
        ggplot2::geom_point(
          ggplot2::aes(color = .data$variant_type_label),
          size = variant_marker_size,
          shape = 17
        )
      p_variant <- apply_variant_marker_color_scale(
        p_variant,
        variant_types = vars$variant_type_label,
        variant_palette = variant_palette,
        variant_colors = variant_colors,
        alpha = variant_alpha,
        name = "Variant type",
        drop = FALSE
      )
    }
    p_variant <- p_variant +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$variant_label),
        angle = 45,
        hjust = 0,
        vjust = 0.5,
        size = text_size,
        color = "black"
      ) +
      ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = x_limits, expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0.7, 1.45), expand = c(0, 0)) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_void() +
      ggplot2::theme(
        text = ggplot2::element_text(color = "black"),
        legend.position = "none",
        plot.margin = ggplot2::margin(8, 12, 0, 8)
      )
    p <- patchwork::wrap_plots(
      p_variant,
      p_connector,
      p_table,
      ncol = 1,
      heights = c(0.75, connector_height, table_height)
    )
  }

  attr(p, "plot_data") <- table_long[]
  attr(p, "variant_data") <- vars[]
  attr(p, "gene_data") <- gene_data
  p
}

parse_hap_numeric_id <- function(hap_id) {
  x <- as.character(hap_id)
  num <- suppressWarnings(as.integer(sub("^.*?([0-9]+)$", "\\1", x)))
  num[is.na(num)] <- .Machine$integer.max
  num
}

make_hap_table_mapper <- function(vars) {
  pos <- sort(unique(as.numeric(vars[["pos"]])))
  idx <- seq_along(pos)
  function(x) {
    x <- as.numeric(x)
    idx[match(x, pos)]
  }
}

make_hap_gene_mapper <- function(vars, region) {
  pos <- sort(unique(as.numeric(vars[["pos"]])))
  n <- length(pos)
  start <- suppressWarnings(as.numeric(region$start %||% min(pos, na.rm = TRUE)))
  end <- suppressWarnings(as.numeric(region$end %||% max(pos, na.rm = TRUE)))
  if (is.na(start)) start <- min(pos, na.rm = TRUE)
  if (is.na(end)) end <- max(pos, na.rm = TRUE)
  if (start > min(pos, na.rm = TRUE)) start <- min(pos, na.rm = TRUE)
  if (end < max(pos, na.rm = TRUE)) end <- max(pos, na.rm = TRUE)
  if (isTRUE(all.equal(start, end))) {
    start <- start - 1
    end <- end + 1
  }

  function(x) {
    x <- as.numeric(x)
    0.5 + (x - start) / (end - start) * n
  }
}

prepare_hap_gene_track <- function(annotation, hap, vars, mapper) {
  region <- hap$region
  chrom <- as.character(region$chrom)[1L]
  query_start <- as.integer(region$start)[1L]
  query_end <- as.integer(region$end)[1L]

  feature <- tryCatch(
    retrieve_feature(annotation, chrom = chrom, start = query_start, end = query_end, as = "Feature"),
    error = function(e) NULL
  )
  if (is.null(feature) || is.null(feature$transcripts) || is.null(feature$exons)) {
    return(NULL)
  }

  tx <- data.table::as.data.table(feature$transcripts)
  ex <- data.table::as.data.table(feature$exons)
  if (nrow(tx) == 0L || nrow(ex) == 0L) {
    return(NULL)
  }

  tx <- tx[as.character(tx[["chrom"]]) == chrom]
  ex <- ex[as.character(ex[["chrom"]]) == chrom]
  tx <- tx[as.integer(tx[["tx_start"]]) <= query_end & as.integer(tx[["tx_end"]]) >= query_start]
  ex <- ex[as.integer(ex[["exon_start"]]) <= query_end & as.integer(ex[["exon_end"]]) >= query_start]
  if (nrow(tx) == 0L || nrow(ex) == 0L) {
    return(NULL)
  }

  if (!is.null(region$locator) && region$locator == "gene" && "gene_id" %in% names(tx)) {
    region_id <- as.character(region$id)[1L]
    tx2 <- tx[as.character(tx[["gene_id"]]) == region_id]
    if (nrow(tx2) > 0L) tx <- tx2
  }
  if (!is.null(region$locator) && region$locator == "transcript" && "transcript_id" %in% names(tx)) {
    region_id <- as.character(region$id)[1L]
    tx2 <- tx[as.character(tx[["transcript_id"]]) == region_id]
    if (nrow(tx2) > 0L) tx <- tx2
  }

  tx_ids <- as.character(tx[["transcript_id"]])
  ex <- ex[as.character(ex[["transcript_id"]]) %in% tx_ids]
  if (nrow(ex) == 0L) {
    return(NULL)
  }

  tx[, ".tx_span" := pmax(1L, as.integer(tx[["tx_end"]]) - as.integer(tx[["tx_start"]]) + 1L)]
  data.table::setorderv(tx, c("gene_id", ".tx_span", "tx_start", "tx_end", "transcript_id"))
  tx[, "track_y" := rev(seq_len(.N))]
  tx[, "tx_xstart" := mapper(pmax(as.integer(tx_start), query_start))]
  tx[, "tx_xend" := mapper(pmin(as.integer(tx_end), query_end))]

  y_map <- tx[, .(transcript_id, track_y)]
  seg <- make_feature_segments(tx, ex)
  if (nrow(seg) == 0L) {
    return(NULL)
  }
  seg <- merge(seg, y_map, by = "transcript_id", all.x = TRUE, sort = FALSE)
  seg <- seg[!is.na(track_y)]
  seg[, "plot_start" := mapper(pmax(query_start, as.integer(get("start"))))]
  seg[, "plot_end" := mapper(pmin(query_end, as.integer(get("end"))))]
  seg <- seg[!is.na(plot_start) & !is.na(plot_end) & plot_start <= plot_end]

  axis_breaks <- make_hap_gene_position_breaks(
    chrom = chrom,
    start = query_start,
    end = query_end,
    mapper = mapper,
    n = 5L
  )

  list(
    transcripts = tx[],
    segments = seg[],
    region = list(chrom = chrom, start = query_start, end = query_end),
    axis_breaks = axis_breaks
  )
}


make_hap_gene_position_breaks <- function(chrom,
                                         start,
                                         end,
                                         mapper,
                                         n = 5L) {
  start <- suppressWarnings(as.numeric(start)[1L])
  end <- suppressWarnings(as.numeric(end)[1L])
  n <- as.integer(n)[1L]
  if (is.na(n) || n < 2L) n <- 5L
  if (is.na(start) || is.na(end)) {
    return(data.table::data.table(x = numeric(), pos = numeric(), label = character()))
  }
  if (start > end) {
    tmp <- start
    start <- end
    end <- tmp
  }
  if (isTRUE(all.equal(start, end))) {
    positions <- start
  } else {
    positions <- pretty(c(start, end), n = n)
    positions <- positions[positions >= start & positions <= end]
    if (!start %in% positions) positions <- c(start, positions)
    if (!end %in% positions) positions <- c(positions, end)
    positions <- sort(unique(round(positions)))
  }
  x <- mapper(positions)
  keep <- !is.na(x) & is.finite(x)
  positions <- positions[keep]
  x <- x[keep]
  if (length(positions) == 0L) {
    return(data.table::data.table(x = numeric(), pos = numeric(), label = character()))
  }
  labels <- format(as.integer(positions), big.mark = ",", scientific = FALSE, trim = TRUE)
  data.table::data.table(x = x, pos = positions, label = labels)
}

draw_hap_gene_track <- function(gene_data,
                                vars,
                                x_limits,
                                x_breaks,
                                x_labels,
                                text_size = 14,
                                gene_text_size = NULL,
                                exon_height = 0.22,
                                cds_height = 0.44,
                                gene_palette = "Paired",
                                gene_colors = NULL,
                                gene_border_color = NA,
                                variant_palette = "Paired",
                                variant_colors = NULL,
                                variant_alpha = 0.6,
                                show_variant_marker = TRUE,
                                variant_marker_size = 2.8,
                                gene_track_legend_position = c("right", "top", "none"),
                                direction_mode = c("transcript", "gene", "end", "none"),
                                show_gene_pos_axis = TRUE,
                                gene_pos_axis_n = 5L,
                                gene_pos_axis_label = "Genomic position (bp)",
                                gene_pos_x_angle = 0) {
  tx <- data.table::copy(gene_data$transcripts)
  seg <- data.table::copy(gene_data$segments)
  vars <- data.table::copy(vars)
  max_y <- max(tx$track_y, na.rm = TRUE)
  if (is.null(gene_text_size)) {
    gene_text_size <- text_size
  }
  gene_text_size <- as.numeric(gene_text_size)[1L]
  if (is.na(gene_text_size) || gene_text_size <= 0) {
    gene_text_size <- text_size
  }
  gene_pos_x_angle <- as.numeric(gene_pos_x_angle)[1L]
  if (is.na(gene_pos_x_angle)) gene_pos_x_angle <- 0
  gene_pos_x_angle <- max(0, min(180, gene_pos_x_angle))
  gene_track_legend_position <- match.arg(gene_track_legend_position)
  direction_mode <- match.arg(direction_mode)
  show_variant_marker <- isTRUE(show_variant_marker)
  variant_marker_size <- as.numeric(variant_marker_size)[1L]
  if (is.na(variant_marker_size)) variant_marker_size <- 2.8
  if (variant_marker_size <= 0) {
    show_variant_marker <- FALSE
    variant_marker_size <- 0
  }
  axis_breaks <- gene_data$axis_breaks
  if (isTRUE(show_gene_pos_axis) && !is.null(gene_data$region)) {
    axis_breaks <- make_hap_gene_position_breaks(
      chrom = gene_data$region$chrom,
      start = gene_data$region$start,
      end = gene_data$region$end,
      mapper = function(x) {
        start <- as.numeric(gene_data$region$start)
        end <- as.numeric(gene_data$region$end)
        n_span <- diff(x_limits)
        if (isTRUE(all.equal(start, end))) {
          return(rep(mean(x_limits), length(x)))
        }
        x_limits[1L] + (as.numeric(x) - start) / (end - start) * n_span
      },
      n = gene_pos_axis_n
    )
  }

  seg[, "feature" := normalize_gene_model_feature(seg[["feature"]])]
  seg[, "feature_width" := ifelse(as.character(feature) == "CDS", cds_height, exon_height)]
  seg[, "ymin" := as.numeric(track_y) - as.numeric(feature_width) / 2]
  seg[, "ymax" := as.numeric(track_y) + as.numeric(feature_width) / 2]

  direction_dt <- make_hap_direction_arrows(
    tx = tx,
    direction_mode = direction_mode,
    model_height = max(exon_height, cds_height)
  )

  vars <- assign_variant_gene_track_y(
    vars = vars,
    tx = tx,
    cds_height = cds_height
  )
  y_upper <- max(seg$ymax, na.rm = TRUE) + 0.35
  if (nrow(direction_dt) > 0L) {
    y_upper <- max(y_upper, max(direction_dt$y, na.rm = TRUE) + 0.15)
  }
  y_lower <- min(c(0.5, vars$marker_y - 0.18), na.rm = TRUE)
  vars[, "connector_top_y" := y_lower]

  x_axis_breaks <- NULL
  x_axis_labels <- NULL
  x_axis_title <- NULL
  x_axis_text <- ggplot2::element_blank()
  x_axis_ticks <- ggplot2::element_blank()
  if (isTRUE(show_gene_pos_axis) && !is.null(axis_breaks) && nrow(axis_breaks) > 0L) {
    x_axis_breaks <- axis_breaks$x
    x_axis_labels <- axis_breaks$label
    x_axis_title <- gene_pos_axis_label
    if (is.null(x_axis_title) || is.na(x_axis_title) || x_axis_title == "") {
      chr_label <- as.character(gene_data$region$chrom %||% NA_character_)
      if (!is.na(chr_label) && nzchar(chr_label)) {
        x_axis_title <- paste0("Chromosome ", chr_label, " position (bp)")
      } else {
        x_axis_title <- "Genomic position (bp)"
      }
    }
    axis_hjust <- if (gene_pos_x_angle == 0) 0.5 else if (gene_pos_x_angle == 90) 1 else 0
    axis_vjust <- if (gene_pos_x_angle == 0) 0 else 0.5
    x_axis_text <- ggplot2::element_text(
      size = gene_text_size,
      angle = gene_pos_x_angle,
      hjust = axis_hjust,
      vjust = axis_vjust,
      color = "black"
    )
    x_axis_ticks <- ggplot2::element_line(color = "black", linewidth = 0.25)
  }

  fill_scale <- make_hap_variant_fill_scale(
    gene_features = unique(as.character(seg$feature)),
    variant_types = unique(as.character(vars$variant_type_label)),
    gene_palette = gene_palette,
    gene_colors = gene_colors,
    variant_palette = variant_palette,
    variant_colors = variant_colors,
    variant_alpha = variant_alpha
  )
  opaque_variant_colors <- make_variant_marker_fill_colors(
    variant_palette = variant_palette,
    variant_colors = variant_colors,
    variant_types = fill_scale$variant_breaks,
    alpha = 1
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = tx,
      ggplot2::aes(x = .data$tx_xstart, xend = .data$tx_xend, y = .data$track_y, yend = .data$track_y),
      linewidth = 0.35,
      color = "black"
    ) +
    ggplot2::geom_segment(
      data = direction_dt,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
      arrow = grid::arrow(length = grid::unit(0.08, "inches")),
      linewidth = 0.25,
      color = "black",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_rect(
      data = seg,
      ggplot2::aes(
        xmin = .data$plot_start,
        xmax = .data$plot_end,
        ymin = .data$ymin,
        ymax = .data$ymax,
        fill = .data$feature
      ),
      color = gene_border_color,
      linewidth = 0.2
    ) +
    ggplot2::geom_segment(
      data = vars,
      ggplot2::aes(x = .data$gene_x, xend = .data$gene_x, y = .data$marker_y, yend = .data$connector_top_y),
      linewidth = 0.35,
      color = "grey35",
      lineend = "butt"
    )
  if (isTRUE(show_variant_marker)) {
    p <- p +
      ggplot2::geom_point(
        data = vars,
        ggplot2::aes(
          x = .data$gene_x,
          y = .data$marker_y,
          color = .data$variant_type_label
        ),
        size = variant_marker_size,
        shape = 17,
        alpha = variant_alpha,
        show.legend = TRUE,
        key_glyph = "point"
      )
  }
  p <- p +
    ggplot2::scale_fill_manual(
      values = fill_scale$colors,
      guide = "none",
      drop = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = opaque_variant_colors,
      limits = fill_scale$variant_breaks,
      breaks = fill_scale$variant_breaks,
      name = "Variant type",
      drop = FALSE,
      guide = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          shape = 17,
          size = max(3, variant_marker_size),
          alpha = 1,
          stroke = 0.5
        )
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_axis_breaks,
      labels = x_axis_labels,
      limits = x_limits,
      expand = c(0, 0),
      position = "top"
    ) +
    ggplot2::scale_y_continuous(
      breaks = tx$track_y,
      labels = tx$transcript_id,
      limits = c(y_lower, y_upper),
      expand = c(0, 0)
    ) +
    ggplot2::labs(x = x_axis_title, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = x_axis_text,
      axis.ticks.x = x_axis_ticks,
      axis.title.x = ggplot2::element_text(size = gene_text_size, color = "black"),
      text = ggplot2::element_text(color = "black"),
      axis.text.y = ggplot2::element_text(size = gene_text_size, color = "black"),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = ifelse(gene_track_legend_position == "none", "none", gene_track_legend_position),
      legend.direction = ifelse(gene_track_legend_position == "top", "horizontal", "vertical"),
      legend.box = ifelse(gene_track_legend_position == "top", "horizontal", "vertical"),
      legend.title = ggplot2::element_text(size = gene_text_size, color = "black"),
      legend.text = ggplot2::element_text(size = gene_text_size, color = "black"),
      plot.margin = ggplot2::margin(8, 12, 0, 8)
    )

  p
}


make_hap_direction_arrows <- function(tx,
                                      direction_mode = c("transcript", "gene", "end", "none"),
                                      model_height = 0.44,
                                      offset = 0.12) {
  direction_mode <- match.arg(direction_mode)
  tx <- data.table::copy(data.table::as.data.table(tx))
  if (identical(direction_mode, "none") || nrow(tx) == 0L || !"strand" %in% names(tx)) {
    return(data.table::data.table(x = numeric(), xend = numeric(), y = numeric(), strand = character()))
  }

  tx[, "strand" := as.character(strand)]
  tx <- tx[strand %in% c("+", "-")]
  if (nrow(tx) == 0L) {
    return(data.table::data.table(x = numeric(), xend = numeric(), y = numeric(), strand = character()))
  }

  arrow_offset <- as.numeric(model_height)[1L] / 2 + as.numeric(offset)[1L]
  if (!is.finite(arrow_offset)) arrow_offset <- 0.34

  if (identical(direction_mode, "transcript")) {
    out <- tx[, .(
      x = as.numeric(tx_xstart),
      xend = as.numeric(tx_xend),
      y = as.numeric(track_y) + arrow_offset,
      strand = as.character(strand)
    )]
  } else {
    if (!"gene_id" %in% names(tx)) {
      tx[, "gene_id" := as.character(transcript_id)]
    }
    out <- tx[, .(
      x = min(as.numeric(tx_xstart), na.rm = TRUE),
      xend = max(as.numeric(tx_xend), na.rm = TRUE),
      y = mean(as.numeric(track_y), na.rm = TRUE) + arrow_offset,
      strand = as.character(strand[1L])
    ), by = gene_id]
    if (identical(direction_mode, "end")) {
      span <- pmax(0, abs(out$xend - out$x))
      short_len <- pmax(0, pmin(span * 0.15, span))
      plus_idx <- which(out$strand == "+")
      minus_idx <- which(out$strand == "-")
      if (length(plus_idx) > 0L) {
        out$x[plus_idx] <- out$xend[plus_idx] - short_len[plus_idx]
      }
      if (length(minus_idx) > 0L) {
        out$xend[minus_idx] <- out$x[minus_idx] + short_len[minus_idx]
      }
    }
    out[, "gene_id" := NULL]
  }

  minus_idx <- which(out$strand == "-")
  if (length(minus_idx) > 0L) {
    old_x <- out$x[minus_idx]
    out$x[minus_idx] <- out$xend[minus_idx]
    out$xend[minus_idx] <- old_x
  }
  out[is.finite(x) & is.finite(xend) & is.finite(y)][]
}


apply_hap_table_fill_scale <- function(p,
                                       values,
                                       table_palette = "Paired",
                                       table_colors = NULL,
                                       table_alpha = 0.6,
                                       fixed_allele_classes = FALSE) {
  values <- unique(as.character(values))
  values <- values[!is.na(values)]
  if (length(values) == 0L) {
    return(p)
  }

  if (isTRUE(fixed_allele_classes)) {
    levels <- hap_table_fill_levels()
    colors <- make_hap_table_fill_colors(
      table_palette = table_palette,
      table_colors = table_colors,
      table_alpha = table_alpha
    )
    return(
      p + ggplot2::scale_fill_manual(
        values = colors,
        limits = levels,
        breaks = levels,
        drop = FALSE,
        na.value = grDevices::adjustcolor("grey85", alpha.f = table_alpha)
      )
    )
  }

  if (!is.null(table_colors)) {
    table_colors <- as.character(table_colors)
    if (is.null(names(table_colors)) || !any(nzchar(names(table_colors)))) {
      if (length(table_colors) < length(values)) {
        table_colors <- grDevices::colorRampPalette(table_colors)(length(values))
      }
      names(table_colors) <- values
    }
    colors <- table_colors
  } else {
    colors <- normalize_discrete_fill_colors(
      n = length(values),
      color_palette = table_palette,
      fill_colors = NULL
    )
    names(colors) <- values
  }

  colors <- grDevices::adjustcolor(colors, alpha.f = table_alpha)
  p + ggplot2::scale_fill_manual(values = colors, na.value = grDevices::adjustcolor("grey85", alpha.f = table_alpha))
}

hap_table_fill_levels <- function() {
  c("A", "T", "C", "G", "indel")
}

normalize_hap_table_fill_class <- function(genotype_label,
                                           genotype_mode = NA_character_) {
  labels <- as.character(genotype_label)
  if (!identical(as.character(genotype_mode)[1L], "string")) {
    return(labels)
  }

  normalized <- toupper(trimws(labels))
  out <- rep(NA_character_, length(labels))
  is_base <- normalized %in% c("A", "T", "C", "G")
  out[is_base] <- normalized[is_base]
  is_indel <- grepl("^I[0-9]+$", normalized)
  out[is_indel] <- "indel"
  out
}

make_hap_table_fill_colors <- function(table_palette = "Paired",
                                       table_colors = NULL,
                                       table_alpha = 0.6) {
  levels <- hap_table_fill_levels()
  table_alpha <- suppressWarnings(as.numeric(table_alpha)[1L])
  if (!is.finite(table_alpha)) table_alpha <- 0.6
  table_alpha <- max(0, min(1, table_alpha))

  default_colors <- normalize_discrete_fill_colors(
    n = length(levels),
    color_palette = table_palette,
    fill_colors = NULL
  )
  names(default_colors) <- levels

  if (!is.null(table_colors)) {
    supplied <- as.character(table_colors)
    supplied <- supplied[!is.na(supplied) & nzchar(supplied)]
    if (length(supplied) > 0L) {
      supplied_names <- names(supplied)
      if (!is.null(supplied_names) && any(nzchar(supplied_names))) {
        normalized_names <- tolower(trimws(supplied_names))
        normalized_names[normalized_names == "a"] <- "A"
        normalized_names[normalized_names == "t"] <- "T"
        normalized_names[normalized_names == "c"] <- "C"
        normalized_names[normalized_names == "g"] <- "G"
        normalized_names[normalized_names %in% c("ind", "indel")] <- "indel"
        names(supplied) <- normalized_names
        supplied <- supplied[names(supplied) %in% levels]
        supplied <- supplied[!duplicated(names(supplied))]
        default_colors[names(supplied)] <- supplied
      } else {
        supplied <- normalize_discrete_fill_colors(
          n = length(levels),
          color_palette = table_palette,
          fill_colors = supplied
        )
        names(supplied) <- levels
        default_colors <- supplied
      }
    }
  }

  out <- grDevices::adjustcolor(default_colors[levels], alpha.f = table_alpha)
  names(out) <- levels
  out
}

assign_variant_gene_track_y <- function(vars,
                                        tx,
                                        cds_height = 0.44,
                                        marker_offset = 0.12) {
  if (!"variant_type" %in% names(vars)) {
    vars[, "variant_type_label" := "..."]
  } else {
    vars[, "variant_type_label" := normalize_variant_marker_type(variant_type)]
  }

  # Variant markers should be aligned on one baseline. Earlier versions attached
  # each marker to the transcript that covered that variant, which made marker
  # and connector positions jump up and down when a gene had multiple isoforms.
  # Use the bottom-most transcript track as the common anchor for all variants.
  bottom_track_y <- suppressWarnings(min(as.numeric(tx[["track_y"]]), na.rm = TRUE))
  if (!is.finite(bottom_track_y)) bottom_track_y <- 1
  bottom_model_y <- bottom_track_y - cds_height / 2

  vars[, "marker_track_y" := bottom_track_y]
  vars[, "gene_model_top_y" := bottom_track_y + cds_height / 2]
  vars[, "gene_model_bottom_y" := bottom_model_y]
  vars[, "marker_y" := bottom_model_y - marker_offset]
  vars[]
}

format_hap_table_genotype_label <- function(genotype,
                                            ref = NULL,
                                            alt = NULL,
                                            genotype_mode = NA_character_) {
  genotype <- as.character(genotype)
  out <- genotype
  out[is.na(out) | out == "" | out == "."] <- "NA"

  if (!identical(as.character(genotype_mode)[1L], "string")) {
    return(out)
  }

  n <- length(out)
  if (is.null(ref)) ref <- rep(NA_character_, n)
  if (is.null(alt)) alt <- rep(NA_character_, n)
  ref <- rep(as.character(ref), length.out = n)
  alt <- rep(as.character(alt), length.out = n)

  vapply(seq_len(n), function(i) {
    value <- out[i]
    if (is.na(value) || value == "NA" || value == "") {
      return("NA")
    }

    parts <- unlist(strsplit(value, "[|/]", perl = TRUE), use.names = FALSE)
    parts <- parts[!is.na(parts) & nzchar(parts) & parts != "."]
    if (length(parts) == 0L) {
      return("NA")
    }

    compact_parts <- format_hap_allele(parts)
    compact_parts <- compact_parts[!is.na(compact_parts) & nzchar(compact_parts)]
    if (length(compact_parts) == 0L) {
      return("NA")
    }

    ref_compact <- format_hap_allele(ref[i])
    alt_parts <- strsplit(alt[i], ",", fixed = TRUE)[[1L]]
    alt_compact <- format_hap_allele(alt_parts)
    alt_compact <- alt_compact[!is.na(alt_compact) & nzchar(alt_compact)]

    alt_hit <- compact_parts[compact_parts %in% alt_compact]
    if (length(alt_hit) > 0L) {
      return(alt_hit[1L])
    }

    if (!is.na(ref_compact) && ref_compact %in% compact_parts) {
      return(ref_compact)
    }

    compact_parts[1L]
  }, character(1L))
}

make_hap_variant_fill_scale <- function(gene_features,
                                        variant_types,
                                        gene_palette = "Paired",
                                        gene_colors = NULL,
                                        variant_palette = "Paired",
                                        variant_colors = NULL,
                                        variant_alpha = 0.6) {
  gene_features <- unique(normalize_gene_model_feature(gene_features))
  gene_features <- gene_features[!is.na(gene_features) & gene_features != ""]
  variant_types <- unique(normalize_variant_marker_type(variant_types))
  variant_types <- variant_types[!is.na(variant_types) & variant_types != ""]

  gene_cols <- make_gene_model_fill_colors(
    color_palette = gene_palette,
    gene_colors = gene_colors,
    features = gene_features
  )
  gene_cols <- gene_cols[intersect(names(gene_cols), unique(c(gene_model_feature_levels(), gene_features)))]

  var_cols <- make_variant_marker_fill_colors(
    variant_palette = variant_palette,
    variant_colors = variant_colors,
    variant_types = variant_types,
    alpha = variant_alpha
  )
  variant_breaks <- intersect(variant_marker_levels(), variant_types)
  if (length(variant_breaks) == 0L) {
    variant_breaks <- variant_marker_levels()
  }

  list(
    colors = c(gene_cols, var_cols),
    variant_colors = var_cols,
    variant_breaks = variant_breaks
  )
}

#' Plot phenotype values grouped by haplotype
#'
#' @description
#' Merges haplotype assignments and phenotype values, then draws phenotype
#' distributions for each haplotype. Haplotype groups are ordered by sample
#' number from left to right. Fill colors are mapped to the median phenotype
#' value of each haplotype group.
#'
#' @param hap A HapVariant object from `hap_variant()`.
#' @param phenotype A phenotype table returned by `read_pheno()` or a compatible data.frame.
#' @param traits Phenotype trait names. If NULL, all numeric traits are used.
#' @param sample_col Sample column name in phenotype table.
#' @param min_hap_samples Minimum sample number required for a haplotype group.
#' @param plot_type Plot type. One of `violin`, `boxplot`, or `violin_boxplot`.
#' @param test_method Pairwise test method. One of `t.test`, `wilcox.test`, or `ks.test`.
#' @param p_adjust P-value adjustment method passed to `p.adjust()`.
#' @param p_label P-value label style. One of `stars`, `number`, or `both`.
#' @param p_cutoff Significance cutoff used for displaying pairwise comparisons.
#' @param p_value_type Which p-value is used for filtering and labeling. One of `raw` or `adjusted`.
#' @param show_signif_only Logical. Whether to only display significant comparisons.
#' @param show_points Logical. Whether to show sample points.
#' @param x_text_angle Rotation angle for haplotype labels on the x-axis.
#' @param facet_ncol Maximum number of facet columns when multiple traits are plotted.
#' @param strip_label_width Maximum character width for wrapping long facet strip labels.
#' @param strip_text_lineheight Line height for wrapped facet strip labels.
#' @param strip_fill Facet strip background fill color.
#' @param strip_border_color Facet strip border color. Use NULL to remove the border.
#' @param show_outliers Logical. Whether to show boxplot outliers.
#' @param fill_palette RColorBrewer palette name used for median-based haplotype fills.
#' @param fill_colors Optional custom fill colors for haplotypes.
#' @param fill_alpha Alpha value for violin/boxplot fill colors.
#' @param violin_width Violin plot width.
#' @param box_width Boxplot width.
#' @param bracket_step Fraction of y-range used to separate significance brackets.
#' @param bracket_tip_fraction Fraction of bracket vertical spacing used for the short downward bracket tips.
#' @param text_size Text size.
#' @return A list with `figure` and `pvalue` elements. Additional elements include `summary`, `bracket`, and `plot_data`.
#' @examples
#' vcf_file <- system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR")
#' anno_file <- system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR")
#' pheno_file <- system.file("extdata", "gtr_demo_pheno.tsv", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file, mode = "memory", verbose = FALSE)
#' anno <- read_genepred(anno_file, format = "genePredExt", verbose = FALSE)
#' hap <- hap_variant(vcf, annotation = anno, gene_id = "GeneA", genotype_mode = "code", min_variant_number = 1)
#' pheno <- read_pheno(pheno_file, verbose = FALSE)
#' plot_hap_pheno(hap, phenotype = pheno, traits = "seed_weight", min_hap_samples = 3)
#' plot_hap_pheno(hap, phenotype = pheno, traits = "seed_weight", min_hap_samples = 3,
#'                test_method = "wilcox.test", p_label = "number")
#' @export
plot_hap_pheno <- function(hap,
                           phenotype,
                           traits = NULL,
                           sample_col = "sample_id",
                           min_hap_samples = 2L,
                           plot_type = c("violin_boxplot", "violin", "boxplot"),
                           test_method = c("t.test", "wilcox.test", "ks.test"),
                           p_adjust = "BH",
                           p_label = c("stars", "number", "both"),
                           p_cutoff = 0.05,
                           p_value_type = c("raw", "adjusted"),
                           show_signif_only = TRUE,
                           show_points = FALSE,
                           show_outliers = FALSE,
                           fill_palette = "Paired",
                           fill_colors = NULL,
                           fill_alpha = 0.75,
                           violin_width = 0.9,
                           box_width = 0.18,
                           bracket_step = 0.08,
                           bracket_tip_fraction = 0.12,
                           x_text_angle = 90,
                           facet_ncol = 3L,
                           strip_label_width = 24,
                           strip_text_lineheight = 0.9,
                           strip_fill = "white",
                           strip_border_color = NULL,
                           text_size = 14) {
  stop_if_not(inherits(hap, "HapVariant"), "`hap` must be a HapVariant object.")
  plot_type <- match.arg(plot_type)
  test_method <- match.arg(test_method)
  p_label <- match.arg(p_label)
  p_value_type <- match.arg(p_value_type)

  min_hap_samples <- as.integer(min_hap_samples)[1L]
  stop_if_not(!is.na(min_hap_samples) && min_hap_samples >= 1L, "`min_hap_samples` must be a positive integer.")
  p_cutoff <- as.numeric(p_cutoff)[1L]
  if (is.na(p_cutoff) || p_cutoff <= 0 || p_cutoff > 1) {
    p_cutoff <- 0.05
  }
  fill_alpha <- as.numeric(fill_alpha)[1L]
  if (is.na(fill_alpha)) fill_alpha <- 0.75
  fill_alpha <- max(0, min(1, fill_alpha))
  x_text_angle <- as.numeric(x_text_angle)[1L]
  if (is.na(x_text_angle)) x_text_angle <- 90
  facet_ncol <- as.integer(facet_ncol)[1L]
  if (is.na(facet_ncol) || facet_ncol < 1L) facet_ncol <- 3L
  strip_label_width <- as.integer(strip_label_width)[1L]
  if (is.na(strip_label_width) || strip_label_width < 1L) strip_label_width <- 24L
  strip_text_lineheight <- as.numeric(strip_text_lineheight)[1L]
  if (is.na(strip_text_lineheight) || strip_text_lineheight <= 0) strip_text_lineheight <- 0.9
  strip_fill <- as.character(strip_fill)[1L]
  if (is.na(strip_fill) || strip_fill == "") strip_fill <- "white"
  if (is.null(strip_border_color) || length(strip_border_color) == 0L || is.na(strip_border_color[1L]) || strip_border_color[1L] == "") {
    strip_border_color <- NA_character_
  } else {
    strip_border_color <- as.character(strip_border_color)[1L]
  }

  pheno <- data.table::as.data.table(phenotype)
  stop_if_not(sample_col %in% names(pheno), paste0("Sample column was not found: ", sample_col))
  data.table::setnames(pheno, sample_col, "sample_id")

  hap_sample <- data.table::as.data.table(hap$sample_haplotypes)
  stop_if_not(all(c("sample_id", "hap_id") %in% names(hap_sample)), "`hap$sample_haplotypes` must contain sample_id and hap_id columns.")

  dt <- merge(hap_sample[, .(sample_id, hap_id)], pheno, by = "sample_id", all.x = FALSE)
  stop_if_not(nrow(dt) > 0L, "No matched samples were found between haplotypes and phenotype table.")

  hap_count <- dt[, .(sample_n = .N), by = hap_id]
  hap_keep <- hap_count[sample_n >= min_hap_samples, hap_id]
  dt <- dt[hap_id %in% hap_keep]
  hap_count <- hap_count[hap_id %in% hap_keep]
  stop_if_not(length(unique(dt$hap_id)) >= 2L, "At least two haplotypes with enough samples are required.")

  if (is.null(traits)) {
    trait_info <- summary_pheno(dt, sample_col = "sample_id")
    traits <- trait_info[type == "numeric" & !trait %in% c("hap_id", "sample_n"), trait]
  }
  traits <- as.character(traits)
  stop_if_not(length(traits) > 0L, "No numeric phenotype traits were selected.")
  stop_if_not(all(traits %in% names(dt)), "Some traits were not found in phenotype table.")

  long <- data.table::melt(
    dt,
    id.vars = c("sample_id", "hap_id"),
    measure.vars = traits,
    variable.name = "trait",
    value.name = "value",
    variable.factor = FALSE
  )
  long[, "value" := suppressWarnings(as.numeric(value))]
  long <- long[!is.na(value)]
  stop_if_not(nrow(long) > 0L, "No non-missing phenotype values were available.")

  # Haplotype order is controlled by group size from left to right.
  hap_count <- hap_count[order(-sample_n, hap_id)]
  hap_order <- hap_count$hap_id
  hap_axis_labels <- paste0(as.character(hap_count$hap_id), " (", hap_count$sample_n, ")")
  names(hap_axis_labels) <- as.character(hap_count$hap_id)
  long[, "hap_id" := factor(hap_id, levels = hap_order)]

  summary_dt <- long[, .(
    sample_n = .N,
    median_value = stats::median(value, na.rm = TRUE),
    mean_value = mean(value, na.rm = TRUE),
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE)
  ), by = .(trait, hap_id)]

  fill_dt <- make_hap_pheno_fill_table(summary_dt, fill_palette = fill_palette, fill_colors = fill_colors, alpha = fill_alpha)
  long <- merge(long, fill_dt[, .(trait, hap_id, hap_fill)], by = c("trait", "hap_id"), all.x = TRUE, sort = FALSE)

  test_dt <- pairwise_hap_test(long, method = test_method, p_adjust = p_adjust)
  bracket_dt <- make_hap_brackets(
    long = long,
    test_dt = test_dt,
    p_cutoff = p_cutoff,
    p_label = p_label,
    p_value_type = p_value_type,
    show_signif_only = show_signif_only,
    bracket_step = bracket_step,
    bracket_tip_fraction = bracket_tip_fraction
  )

  # Preserve the exact trait order supplied by the user in facet layout.
  plot_long <- data.table::copy(long)
  plot_long[, "trait" := factor(as.character(trait), levels = traits)]
  plot_bracket <- data.table::copy(bracket_dt)
  if (nrow(plot_bracket) > 0L) {
    plot_bracket[, "trait" := factor(as.character(trait), levels = traits)]
  }

  p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = .data$hap_id, y = .data$value, fill = .data$hap_fill))

  if (plot_type %in% c("violin", "violin_boxplot")) {
    p <- p + ggplot2::geom_violin(
      width = violin_width,
      trim = FALSE,
      color = "black",
      linewidth = 0.35,
      alpha = fill_alpha
    )
  }

  if (plot_type %in% c("boxplot", "violin_boxplot")) {
    outlier_shape <- if (isTRUE(show_outliers)) 19 else NA
    p <- p + ggplot2::geom_boxplot(
      width = box_width,
      outlier.shape = outlier_shape,
      color = "black",
      linewidth = 0.35,
      alpha = min(1, fill_alpha + 0.15)
    )
  }

  if (isTRUE(show_points)) {
    p <- p + ggplot2::geom_jitter(
      width = 0.12,
      height = 0,
      size = 1.4,
      alpha = 0.65,
      color = "black"
    )
  }

  if (length(traits) > 1L) {
    p <- p + ggplot2::facet_wrap(
      ggplot2::vars(.data$trait),
      scales = "free_y",
      ncol = min(facet_ncol, length(traits)),
      labeller = ggplot2::as_labeller(function(x) wrap_strip_labels(x, width = strip_label_width))
    )
  }

  if (nrow(bracket_dt) > 0L) {
    p <- p +
      ggplot2::geom_segment(
        data = plot_bracket,
        ggplot2::aes(x = .data$x1, xend = .data$x2, y = .data$y, yend = .data$y),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::geom_segment(
        data = plot_bracket,
        ggplot2::aes(x = .data$x1, xend = .data$x1, y = .data$y, yend = .data$y_tip),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::geom_segment(
        data = plot_bracket,
        ggplot2::aes(x = .data$x2, xend = .data$x2, y = .data$y, yend = .data$y_tip),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::geom_text(
        data = plot_bracket,
        ggplot2::aes(x = .data$x_mid, y = .data$label_y, label = .data$label),
        inherit.aes = FALSE,
        color = "black",
        size = text_size / 3.2,
        vjust = 0
      )
  }

  fill_values <- unique(long[, .(hap_fill)])$hap_fill
  names(fill_values) <- fill_values

  plot_title <- if (length(traits) == 1L) as.character(traits[1L]) else NULL

  p <- p +
    ggplot2::scale_x_discrete(labels = hap_axis_labels) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(x = "Haplotype", y = "Phenotype value", title = plot_title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(size = text_size, angle = x_text_angle, hjust = 1, vjust = 0.5, color = "black"),
      axis.text.y = ggplot2::element_text(size = text_size, color = "black"),
      axis.title = ggplot2::element_text(size = text_size, color = "black"),
      plot.title = ggplot2::element_text(size = text_size, color = "black", hjust = 0.5),
      strip.text = ggplot2::element_text(size = text_size, color = "black", lineheight = strip_text_lineheight, margin = ggplot2::margin(3, 3, 3, 3)),
      strip.background = ggplot2::element_rect(fill = strip_fill, color = strip_border_color),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "none"
    )

  attr(p, "plot_data") <- long[]
  attr(p, "summary_table") <- summary_dt[]
  attr(p, "test_table") <- test_dt[]
  attr(p, "bracket_table") <- bracket_dt[]

  out <- list(
    figure = p,
    pvalue = test_dt[],
    summary = summary_dt[],
    bracket = bracket_dt[],
    plot_data = long[]
  )
  class(out) <- c("GeneTrackRPhenoPlot", "list")
  out
}



wrap_strip_labels <- function(labels, width = 24) {
  width <- as.integer(width)[1L]
  if (is.na(width) || width < 1L) width <- 24L
  vapply(as.character(labels), function(label) {
    if (is.na(label) || label == "") {
      return(label)
    }
    label <- gsub("[_\\.]+", " ", label)
    label <- gsub("\\s+", " ", trimws(label))
    pieces <- unlist(strsplit(label, " ", fixed = TRUE), use.names = FALSE)
    pieces <- pieces[nzchar(pieces)]
    if (length(pieces) == 0L) {
      return(label)
    }
    chunks <- character()
    current <- ""
    for (piece in pieces) {
      if (nchar(piece, type = "width") > width) {
        starts <- seq(1L, nchar(piece), by = width)
        hard <- substring(piece, starts, pmin(starts + width - 1L, nchar(piece)))
        if (nzchar(current)) {
          chunks <- c(chunks, current)
          current <- ""
        }
        chunks <- c(chunks, hard)
      } else if (!nzchar(current)) {
        current <- piece
      } else if (nchar(paste(current, piece), type = "width") <= width) {
        current <- paste(current, piece)
      } else {
        chunks <- c(chunks, current)
        current <- piece
      }
    }
    if (nzchar(current)) {
      chunks <- c(chunks, current)
    }
    paste(chunks, collapse = "\n")
  }, character(1L), USE.NAMES = FALSE)
}

make_hap_pheno_fill_table <- function(summary_dt,
                                      fill_palette = "Paired",
                                      fill_colors = NULL,
                                      alpha = 0.75) {
  dt <- data.table::copy(summary_dt)
  dt[, "hap_fill" := NA_character_]

  for (trait_i in unique(as.character(dt$trait))) {
    idx <- which(as.character(dt$trait) == trait_i)
    med <- dt$median_value[idx]
    n <- length(idx)
    if (!is.null(fill_colors)) {
      cols <- as.character(fill_colors)
      if (!is.null(names(cols)) && all(as.character(dt$hap_id[idx]) %in% names(cols))) {
        cols <- cols[as.character(dt$hap_id[idx])]
      } else {
        if (length(cols) < n) cols <- grDevices::colorRampPalette(cols)(n)
        cols <- cols[seq_len(n)]
      }
    } else {
      cols <- normalize_discrete_fill_colors(n = n, color_palette = fill_palette, fill_colors = NULL)
    }
    ord <- order(med, na.last = TRUE)
    assigned <- rep(NA_character_, n)
    assigned[ord] <- cols[seq_len(n)]
    dt$hap_fill[idx] <- grDevices::adjustcolor(assigned, alpha.f = alpha)
  }

  dt[]
}

pairwise_hap_test <- function(long,
                              method = "t.test",
                              p_adjust = "BH") {
  out <- long[, {
    haps <- as.character(stats::na.omit(unique(hap_id)))
    if (length(haps) < 2L) {
      return(data.table::data.table(group1 = character(), group2 = character(), p_value = numeric()))
    }
    pairs <- utils::combn(haps, 2L, simplify = FALSE)
    res <- lapply(pairs, function(pair) {
      x <- value[as.character(hap_id) == pair[1L]]
      y <- value[as.character(hap_id) == pair[2L]]
      pval <- run_hap_pair_test(x, y, method = method)
      data.table::data.table(group1 = pair[1L], group2 = pair[2L], p_value = pval)
    })
    data.table::rbindlist(res)
  }, by = trait]
  out[, "p_adj" := stats::p.adjust(p_value, method = p_adjust), by = trait]
  out[, "method" := method]
  out[]
}

run_hap_pair_test <- function(x, y, method = "t.test") {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) < 2L || length(y) < 2L) {
    return(NA_real_)
  }
  tryCatch({
    if (identical(method, "t.test")) {
      stats::t.test(x, y)$p.value
    } else if (identical(method, "wilcox.test")) {
      stats::wilcox.test(x, y)$p.value
    } else if (identical(method, "ks.test")) {
      stats::ks.test(x, y)$p.value
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
}

make_hap_brackets <- function(long,
                              test_dt,
                              p_cutoff = 0.05,
                              p_label = "stars",
                              p_value_type = "adjusted",
                              show_signif_only = TRUE,
                              bracket_step = 0.08,
                              bracket_tip_fraction = 0.12) {
  if (nrow(test_dt) == 0L) return(data.table::data.table())
  dt <- data.table::copy(test_dt)
  dt[, "p_display" := if (identical(p_value_type, "raw")) p_value else p_adj]
  if (isTRUE(show_signif_only)) {
    dt <- dt[!is.na(p_display) & p_display <= p_cutoff]
  } else {
    dt <- dt[!is.na(p_display)]
  }
  if (nrow(dt) == 0L) return(data.table::data.table())

  level_dt <- unique(long[, .(hap_id)])
  level_dt[, "x" := as.integer(hap_id)]
  level_dt[, "hap_id_chr" := as.character(hap_id)]

  y_info <- long[, .(
    y_max = max(value, na.rm = TRUE),
    y_min = min(value, na.rm = TRUE)
  ), by = trait]
  y_info[, "y_range" := y_max - y_min]
  y_info[y_range == 0 | is.na(y_range), "y_range" := abs(y_max)]
  y_info[y_range == 0 | is.na(y_range), "y_range" := 1]
  dt <- merge(dt, y_info, by = "trait", all.x = TRUE, sort = FALSE)

  dt[, "x1" := level_dt$x[match(group1, level_dt$hap_id_chr)]]
  dt[, "x2" := level_dt$x[match(group2, level_dt$hap_id_chr)]]
  dt <- dt[!is.na(x1) & !is.na(x2)]
  if (nrow(dt) == 0L) return(data.table::data.table())

  data.table::setorder(dt, trait, p_display, group1, group2)
  dt[, "bracket_index" := seq_len(.N), by = trait]
  step <- as.numeric(bracket_step)[1L]
  if (is.na(step) || step <= 0) step <- 0.08
  tip_fraction <- as.numeric(bracket_tip_fraction)[1L]
  if (is.na(tip_fraction) || tip_fraction < 0) tip_fraction <- 0.12
  tip_fraction <- min(tip_fraction, 1)
  dt[, "y" := y_max + y_range * step * bracket_index]
  dt[, "y_tip" := y - y_range * step * tip_fraction]
  dt[, "label_y" := y + y_range * step * 0.05]
  dt[, "x_mid" := (x1 + x2) / 2]
  dt[, "label" := format_hap_p_label(p_display, style = p_label)]
  dt[]
}

format_hap_p_label <- function(p, style = "stars") {
  stars <- ifelse(
    is.na(p), "NA",
    ifelse(p <= 0.001, "***", ifelse(p <= 0.01, "**", ifelse(p <= 0.05, "*", "ns")))
  )
  num <- format_p_value(p)
  if (identical(style, "number")) {
    num
  } else if (identical(style, "both")) {
    paste0(stars, "\n", num)
  } else {
    stars
  }
}

format_p_value <- function(p) {
  out <- ifelse(
    is.na(p),
    "NA",
    ifelse(
      p < 0.001,
      formatC(p, format = "e", digits = 2),
      sprintf("%.3g", p)
    )
  )
  out
}
