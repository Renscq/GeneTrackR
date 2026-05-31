# Author: Rensc
# Date: 2026-05-29
# Version: 0.1.0
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
#' @param show_reference_row Logical. Whether to add the first table row showing REF/ALT allele strings for each variant.
#' @param variant_label Column used for variant labels. One of `variant_id`, `pos`, or an existing column in `hap$variants`.
#' @param text_size Base text size.
#' @param table_text_size Haplotype table text size.
#' @param table_x_angle Angle of haplotype table x-axis labels. Default is 90.
#' @param genotype_text_size Genotype cell text size. If NULL, `table_text_size` is used.
#' @param gene_track_height Relative height of the gene track panel.
#' @param connector_height Relative height of the connector panel.
#' @param table_height Relative height of the haplotype table panel.
#' @param exon_height Height of exon boxes in the compact gene track.
#' @param cds_height Height of CDS boxes in the compact gene track.
#' @param color_palette RColorBrewer palette name used for gene model feature fills.
#' @param fill_colors Optional custom fill colors for gene model features.
#' @param border_color Border color for gene model feature boxes.
#' @param table_fill_palette RColorBrewer palette name used for haplotype table fills.
#' @param table_fill_colors Optional custom fill colors for haplotype table genotypes.
#' @param table_fill_alpha Alpha value for haplotype table fill colors.
#' @param reference_fill Background fill color for the REF/ALT reference row.
#' @param variant_palette RColorBrewer palette name used for variant-type marker fills.
#' @param variant_colors Optional custom fill colors for variant-type markers.
#' @return A patchwork object with attributes `plot_data`, `variant_data`, and `gene_data`.
#' @examples
#' vcf_file <- system.file("extdata", "example_haplotype.vcf", package = "GeneTrackR")
#' anno_file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file)
#' anno <- read_genepred(anno_file, format = "genePredExt", verbose = FALSE)
#' hap <- hap_variant(vcf, annotation = anno, gene_id = "GeneA", genotype_mode = "string")
#' plot_hap_variant(hap, annotation = anno, min_hap_samples = 1)
#' plot_hap_variant(
#'   hap,
#'   annotation = anno,
#'   min_hap_samples = 1,
#'   table_x_angle = 90,
#'   table_fill_palette = "RdBu"
#' )
#' @export
plot_hap_variant <- function(hap,
                             annotation = NULL,
                             show_gene_model = TRUE,
                             min_hap_samples = 5L,
                             show_reference_row = TRUE,
                             variant_label = c("variant_id", "pos"),
                             text_size = 14,
                             table_text_size = 3.2,
                             table_x_angle = 90,
                             genotype_text_size = NULL,
                             gene_track_height = 1.25,
                             connector_height = 0.35,
                             table_height = NULL,
                             exon_height = 0.22,
                             cds_height = 0.44,
                             color_palette = "Paired",
                             fill_colors = NULL,
                             border_color = NA,
                             table_fill_palette = "RdBu",
                             table_fill_colors = NULL,
                             table_fill_alpha = 0.6,
                             reference_fill = "white",
                             variant_palette = "Set2",
                             variant_colors = NULL) {
  stop_if_not(inherits(hap, "HapVariant"), "`hap` must be a HapVariant object.")

  vars <- data.table::as.data.table(hap$variants)
  haps <- data.table::as.data.table(hap$haplotypes)
  stop_if_not(nrow(vars) > 0L, "No variants are available in `hap$variants`.")
  stop_if_not(nrow(haps) > 0L, "No haplotypes are available in `hap$haplotypes`.")

  min_hap_samples <- as.integer(min_hap_samples)[1L]
  stop_if_not(!is.na(min_hap_samples) && min_hap_samples >= 1L, "`min_hap_samples` must be a positive integer.")
  show_reference_row <- isTRUE(show_reference_row)
  border_color <- normalize_border_color(border_color)
  table_fill_alpha <- as.numeric(table_fill_alpha)[1L]
  if (is.na(table_fill_alpha)) table_fill_alpha <- 0.6
  table_fill_alpha <- max(0, min(1, table_fill_alpha))
  if (is.null(genotype_text_size)) {
    genotype_text_size <- table_text_size
  }
  genotype_text_size <- as.numeric(genotype_text_size)[1L]
  if (is.na(genotype_text_size) || genotype_text_size <= 0) {
    genotype_text_size <- table_text_size
  }
  table_x_angle <- as.numeric(table_x_angle)[1L]
  if (is.na(table_x_angle)) table_x_angle <- 90
  table_x_angle <- max(0, min(180, table_x_angle))
  reference_fill <- as.character(reference_fill)[1L]
  if (is.na(reference_fill) || reference_fill == "") {
    reference_fill <- "white"
  }

  haps <- haps[sample_n >= min_hap_samples]
  stop_if_not(nrow(haps) > 0L, "No haplotypes remained after `min_hap_samples` filtering.")

  variant_ids <- vars[["variant_id"]]
  variant_ids <- variant_ids[variant_ids %in% names(haps)]
  stop_if_not(length(variant_ids) > 0L, "No variant genotype columns were found in `hap$haplotypes`.")

  vars <- vars[match(variant_ids, variant_id)]
  vars[, "variant_index" := seq_len(.N)]

  variant_label <- match.arg(variant_label)
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

  hap_order <- haps[order(-sample_n, hap_id), hap_id]
  table_long[, "row_label" := paste0(hap_id, " (n=", sample_n, ")")]
  hap_label_levels <- haps[match(hap_order, hap_id), paste0(hap_id, " (n=", sample_n, ")")]
  row_levels <- rev(hap_label_levels)

  table_long[, "genotype_label" := format_hap_table_genotype_label(
    genotype = genotype,
    ref = ref,
    alt = alt,
    genotype_mode = hap$meta$genotype_mode %||% NA_character_
  )]
  table_long[is.na(genotype_label), "genotype_label" := "NA"]

  if (show_reference_row) {
    ref_row_label <- "REF/ALT"
    ref_alt <- vars[, .(
      variant_id,
      variant_index,
      variant_label,
      row_label = ref_row_label,
      genotype = format_hap_ref_alt(ref, alt),
      genotype_label = format_hap_ref_alt(ref, alt)
    )]
    for (nm in setdiff(names(table_long), names(ref_alt))) {
      ref_alt[, (nm) := NA]
    }
    data.table::setcolorder(ref_alt, names(table_long))
    table_long <- data.table::rbindlist(list(table_long, ref_alt), fill = TRUE)
    row_levels <- c(row_levels, ref_row_label)
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
        exon_height = exon_height,
        cds_height = cds_height,
        color_palette = color_palette,
        fill_colors = fill_colors,
        border_color = border_color,
        variant_palette = variant_palette,
        variant_colors = variant_colors
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
    table_reference <- table_long[row_label == "REF/ALT"]
    table_genotype <- table_long[row_label != "REF/ALT"]
  }

  p_table <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = table_genotype,
      ggplot2::aes(x = .data$variant_index, y = .data$hap_y, fill = .data$genotype_label),
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
      data = table_long,
      ggplot2::aes(x = .data$variant_index, y = .data$hap_y, label = .data$genotype_label),
      size = genotype_text_size,
      color = "black"
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
      axis.text.x = ggplot2::element_text(size = text_size * 0.68, angle = table_x_angle, hjust = 1, vjust = 0.5, color = "black"),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = text_size * 0.78, color = "black"),
      axis.title.y = ggplot2::element_text(size = text_size, color = "black"),
      axis.title.x = ggplot2::element_blank(),
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 12, 8, 8)
    )

  p_table <- apply_hap_table_fill_scale(
    p_table,
    values = table_genotype$genotype_label,
    table_fill_palette = table_fill_palette,
    table_fill_colors = table_fill_colors,
    table_fill_alpha = table_fill_alpha
  )

  if (is.null(table_height)) {
    table_height <- max(1.2, n_rows * 0.38)
  }

  if (!is.null(p_gene)) {
    p <- p_gene / p_connector / p_table +
      patchwork::plot_layout(heights = c(gene_track_height, connector_height, table_height))
  } else {
    p_variant <- ggplot2::ggplot(vars, ggplot2::aes(x = .data$gene_x, y = 1)) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$variant_label),
        angle = 45,
        hjust = 0,
        vjust = 0.5,
        size = genotype_text_size,
        color = "black"
      ) +
      ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = x_limits, expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0.7, 1.45), expand = c(0, 0)) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_void() +
      ggplot2::theme(
        text = ggplot2::element_text(color = "black"),
        plot.margin = ggplot2::margin(8, 12, 0, 8)
      )
    p <- p_variant / p_connector / p_table +
      patchwork::plot_layout(heights = c(0.75, connector_height, table_height))
  }

  attr(p, "plot_data") <- table_long[]
  attr(p, "variant_data") <- vars[]
  attr(p, "gene_data") <- gene_data
  p
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

  data.table::setorderv(tx, c("gene_id", "tx_start", "tx_end", "transcript_id"))
  tx[, "track_y" := seq_len(.N)]
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

  list(transcripts = tx[], segments = seg[])
}

draw_hap_gene_track <- function(gene_data,
                                vars,
                                x_limits,
                                x_breaks,
                                x_labels,
                                text_size = 14,
                                exon_height = 0.22,
                                cds_height = 0.44,
                                color_palette = "Paired",
                                fill_colors = NULL,
                                border_color = NA,
                                variant_palette = "Set2",
                                variant_colors = NULL) {
  tx <- data.table::copy(gene_data$transcripts)
  seg <- data.table::copy(gene_data$segments)
  vars <- data.table::copy(vars)
  max_y <- max(tx$track_y, na.rm = TRUE)

  seg[, "feature_width" := ifelse(as.character(feature) == "CDS", cds_height, exon_height)]
  seg[, "ymin" := as.numeric(track_y) - as.numeric(feature_width) / 2]
  seg[, "ymax" := as.numeric(track_y) + as.numeric(feature_width) / 2]

  vars <- assign_variant_gene_track_y(
    vars = vars,
    tx = tx,
    cds_height = cds_height
  )
  y_upper <- max(seg$ymax, na.rm = TRUE) + 0.35
  y_lower <- min(c(0.5, vars$marker_y - 0.18), na.rm = TRUE)
  vars[, "connector_top_y" := y_lower]

  fill_scale <- make_hap_variant_fill_scale(
    gene_features = unique(as.character(seg$feature)),
    variant_types = unique(as.character(vars$variant_type_label)),
    color_palette = color_palette,
    fill_colors = fill_colors,
    variant_palette = variant_palette,
    variant_colors = variant_colors
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = tx,
      ggplot2::aes(x = .data$tx_xstart, xend = .data$tx_xend, y = .data$track_y, yend = .data$track_y),
      linewidth = 0.35,
      color = "black"
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
      color = border_color,
      linewidth = 0.2
    ) +
    ggplot2::geom_segment(
      data = vars,
      ggplot2::aes(x = .data$gene_x, xend = .data$gene_x, y = .data$marker_y, yend = .data$connector_top_y),
      linewidth = 0.35,
      color = "grey35",
      lineend = "butt"
    ) +
    ggplot2::geom_point(
      data = vars,
      ggplot2::aes(x = .data$gene_x, y = .data$marker_y, fill = .data$variant_type_label),
      size = 2.8,
      shape = 24,
      color = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = fill_scale$colors,
      breaks = fill_scale$variant_breaks,
      name = "Variant type",
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
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
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      text = ggplot2::element_text(color = "black"),
      axis.text.y = ggplot2::element_text(size = text_size * 0.7, color = "black"),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = text_size * 0.75, color = "black"),
      legend.text = ggplot2::element_text(size = text_size * 0.68, color = "black"),
      plot.margin = ggplot2::margin(8, 12, 0, 8)
    )

  p
}


apply_hap_table_fill_scale <- function(p,
                                       values,
                                       table_fill_palette = "RdBu",
                                       table_fill_colors = NULL,
                                       table_fill_alpha = 0.6) {
  values <- unique(as.character(values))
  values <- values[!is.na(values)]
  if (length(values) == 0L) {
    return(p)
  }

  if (!is.null(table_fill_colors)) {
    table_fill_colors <- as.character(table_fill_colors)
    if (is.null(names(table_fill_colors)) || !any(nzchar(names(table_fill_colors)))) {
      if (length(table_fill_colors) < length(values)) {
        table_fill_colors <- grDevices::colorRampPalette(table_fill_colors)(length(values))
      }
      names(table_fill_colors) <- values
    }
    colors <- table_fill_colors
  } else {
    colors <- normalize_discrete_fill_colors(
      n = length(values),
      color_palette = table_fill_palette,
      fill_colors = NULL
    )
    names(colors) <- values
  }

  colors <- grDevices::adjustcolor(colors, alpha.f = table_fill_alpha)
  p + ggplot2::scale_fill_manual(values = colors, na.value = grDevices::adjustcolor("grey85", alpha.f = table_fill_alpha))
}

assign_variant_gene_track_y <- function(vars,
                                        tx,
                                        cds_height = 0.44,
                                        marker_offset = 0.12) {
  if (!"variant_type" %in% names(vars)) {
    vars[, "variant_type_label" := "variant"]
  } else {
    vars[, "variant_type_label" := as.character(variant_type)]
    vars[is.na(variant_type_label) | variant_type_label == "", "variant_type_label" := "variant"]
  }

  vars[, "marker_track_y" := NA_real_]
  vars[, "gene_model_top_y" := NA_real_]
  vars[, "gene_model_bottom_y" := NA_real_]

  tx_start_col <- if ("tx_start" %in% names(tx)) "tx_start" else "start"
  tx_end_col <- if ("tx_end" %in% names(tx)) "tx_end" else "end"

  for (i in seq_len(nrow(vars))) {
    pos_i <- suppressWarnings(as.integer(vars[["pos"]][i]))
    hit <- tx[as.integer(get(tx_start_col)) <= pos_i & as.integer(get(tx_end_col)) >= pos_i]
    if (nrow(hit) == 0L) {
      dist <- pmin(
        abs(pos_i - as.integer(tx[[tx_start_col]])),
        abs(pos_i - as.integer(tx[[tx_end_col]]))
      )
      hit <- tx[which.min(dist)]
    }
    y_i <- as.numeric(hit[["track_y"]][1L])
    data.table::set(vars, i = i, j = "marker_track_y", value = y_i)
    data.table::set(vars, i = i, j = "gene_model_top_y", value = y_i + cds_height / 2)
    data.table::set(vars, i = i, j = "gene_model_bottom_y", value = y_i - cds_height / 2)
  }

  vars[, "marker_y" := gene_model_bottom_y - marker_offset]
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
                                        color_palette = "Paired",
                                        fill_colors = NULL,
                                        variant_palette = "Set2",
                                        variant_colors = NULL) {
  gene_features <- unique(as.character(gene_features))
  gene_features <- gene_features[!is.na(gene_features) & gene_features != ""]
  variant_types <- unique(as.character(variant_types))
  variant_types <- variant_types[!is.na(variant_types) & variant_types != ""]

  gene_cols <- normalize_discrete_fill_colors(
    n = length(gene_features),
    color_palette = color_palette,
    fill_colors = fill_colors
  )
  names(gene_cols) <- gene_features

  if (!is.null(variant_colors)) {
    variant_colors <- as.character(variant_colors)
    if (is.null(names(variant_colors)) || !any(nzchar(names(variant_colors)))) {
      if (length(variant_colors) < length(variant_types)) {
        variant_colors <- grDevices::colorRampPalette(variant_colors)(length(variant_types))
      }
      names(variant_colors) <- variant_types
    }
    var_cols <- variant_colors
  } else {
    var_cols <- normalize_discrete_fill_colors(
      n = length(variant_types),
      color_palette = variant_palette,
      fill_colors = NULL
    )
    names(var_cols) <- variant_types
  }

  list(
    colors = c(gene_cols, var_cols),
    variant_breaks = variant_types
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
#' @param strip_label_width Maximum character width for wrapping long facet strip labels.
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
#' vcf_file <- system.file("extdata", "example_haplotype.vcf", package = "GeneTrackR")
#' anno_file <- system.file("extdata", "example.genePredExt", package = "GeneTrackR")
#' pheno_file <- system.file("extdata", "example_pheno.tsv", package = "GeneTrackR")
#' vcf <- read_vcf(vcf_file)
#' anno <- read_genepred(anno_file, format = "genePredExt", verbose = FALSE)
#' hap <- hap_variant(vcf, annotation = anno, gene_id = "GeneA", genotype_mode = "code")
#' pheno <- read_pheno(pheno_file)
#' plot_hap_pheno(hap, phenotype = pheno, traits = "plant_height", min_hap_samples = 1)
#' plot_hap_pheno(hap, phenotype = pheno, traits = "plant_height", min_hap_samples = 1,
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
                           fill_palette = "RdBu",
                           fill_colors = NULL,
                           fill_alpha = 0.75,
                           violin_width = 0.9,
                           box_width = 0.18,
                           bracket_step = 0.08,
                           bracket_tip_fraction = 0.12,
                           x_text_angle = 90,
                           strip_label_width = 24,
                           strip_fill = "white",
                           strip_border_color = NULL,
                           strip_text_lineheight = 0.9,
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


  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$hap_id, y = .data$value, fill = .data$hap_fill))

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
      labeller = ggplot2::as_labeller(function(x) wrap_strip_labels(x, width = strip_label_width))
    )
  }

  if (nrow(bracket_dt) > 0L) {
    p <- p +
      ggplot2::geom_segment(
        data = bracket_dt,
        ggplot2::aes(x = .data$x1, xend = .data$x2, y = .data$y, yend = .data$y),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::geom_segment(
        data = bracket_dt,
        ggplot2::aes(x = .data$x1, xend = .data$x1, y = .data$y, yend = .data$y_tip),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::geom_segment(
        data = bracket_dt,
        ggplot2::aes(x = .data$x2, xend = .data$x2, y = .data$y, yend = .data$y_tip),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::geom_text(
        data = bracket_dt,
        ggplot2::aes(x = .data$x_mid, y = .data$label_y, label = .data$label),
        inherit.aes = FALSE,
        color = "black",
        size = text_size / 3.2,
        vjust = 0
      )
  }

  fill_values <- unique(long[, .(hap_fill)])$hap_fill
  names(fill_values) <- fill_values

  p <- p +
    ggplot2::scale_x_discrete(labels = hap_axis_labels) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(x = "Haplotype", y = "Phenotype value") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(size = text_size, angle = x_text_angle, hjust = 1, vjust = 0.5, color = "black"),
      axis.text.y = ggplot2::element_text(size = text_size, color = "black"),
      axis.title = ggplot2::element_text(size = text_size, color = "black"),
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
                                      fill_palette = "RdBu",
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
