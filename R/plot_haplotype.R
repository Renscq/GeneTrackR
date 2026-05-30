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
#' Plot phenotype values grouped by haplotype
#'
#' @description
#' Merges haplotype assignments and phenotype values, then draws boxplots with
#' pairwise t-test summaries.
#'
#' @param hap A HapVariant object from `hap_variant()`.
#' @param phenotype A phenotype table returned by `read_pheno()` or a compatible data.frame.
#' @param traits Phenotype trait names. If NULL, all numeric traits are used.
#' @param sample_col Sample column name in phenotype table.
#' @param min_hap_samples Minimum sample number required for a haplotype group.
#' @param p_adjust P-value adjustment method passed to `p.adjust()`.
#' @param show_points Logical. Whether to show sample points.
#' @param text_size Text size.
#' @return A ggplot object with attributes `plot_data` and `test_table`.
#' @export
plot_hap_pheno <- function(hap,
                           phenotype,
                           traits = NULL,
                           sample_col = "sample_id",
                           min_hap_samples = 2L,
                           p_adjust = "BH",
                           show_points = TRUE,
                           text_size = 14) {
  stop_if_not(inherits(hap, "HapVariant"), "`hap` must be a HapVariant object.")
  pheno <- data.table::as.data.table(phenotype)
  stop_if_not(sample_col %in% names(pheno), paste0("Sample column was not found: ", sample_col))
  data.table::setnames(pheno, sample_col, "sample_id")

  hap_sample <- data.table::as.data.table(hap$sample_haplotypes)
  dt <- merge(hap_sample[, .(sample_id, hap_id)], pheno, by = "sample_id", all.x = FALSE)
  stop_if_not(nrow(dt) > 0L, "No matched samples were found between haplotypes and phenotype table.")

  hap_keep <- dt[, .N, by = hap_id][N >= as.integer(min_hap_samples), hap_id]
  dt <- dt[hap_id %in% hap_keep]
  stop_if_not(length(unique(dt$hap_id)) >= 2L, "At least two haplotypes with enough samples are required.")

  if (is.null(traits)) {
    trait_info <- summary_pheno(dt, sample_col = "sample_id")
    traits <- trait_info[type == "numeric" & trait != "hap_id", trait]
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
  long <- long[!is.na(value)]
  stop_if_not(nrow(long) > 0L, "No non-missing phenotype values were available.")

  test_dt <- pairwise_hap_ttest(long, p_adjust = p_adjust)
  label_dt <- make_hap_test_labels(long, test_dt)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$hap_id, y = .data$value, fill = .data$hap_id)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    ggplot2::facet_wrap(ggplot2::vars(.data$trait), scales = "free_y") +
    ggplot2::labs(x = "Haplotype", y = "Phenotype value") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(size = text_size, angle = 45, hjust = 1, color = "black"),
      axis.text.y = ggplot2::element_text(size = text_size, color = "black"),
      axis.title = ggplot2::element_text(size = text_size, color = "black"),
      strip.text = ggplot2::element_text(size = text_size, color = "black"),
      legend.position = "none"
    )

  if (isTRUE(show_points)) {
    p <- p + ggplot2::geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.75)
  }

  if (nrow(label_dt) > 0L) {
    p <- p + ggplot2::geom_text(
      data = label_dt,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      inherit.aes = FALSE,
      size = text_size / 3.2,
      hjust = 0.5,
      vjust = 1
    )
  }

  attr(p, "plot_data") <- long[]
  attr(p, "test_table") <- test_dt[]
  p
}

pairwise_hap_ttest <- function(long, p_adjust = "BH") {
  out <- long[, {
    haps <- unique(hap_id)
    if (length(haps) < 2L) {
      return(data.table::data.table(group1 = character(), group2 = character(), p_value = numeric()))
    }
    pairs <- utils::combn(haps, 2L, simplify = FALSE)
    res <- lapply(pairs, function(pair) {
      x <- value[hap_id == pair[1L]]
      y <- value[hap_id == pair[2L]]
      if (length(x) < 2L || length(y) < 2L) {
        pval <- NA_real_
      } else {
        pval <- tryCatch(stats::t.test(x, y)$p.value, error = function(e) NA_real_)
      }
      data.table::data.table(group1 = pair[1L], group2 = pair[2L], p_value = pval)
    })
    data.table::rbindlist(res)
  }, by = trait]
  out[, "p_adj" := stats::p.adjust(p_value, method = p_adjust), by = trait]
  out[, "label" := paste0(group1, " vs ", group2, ": p=", format_p_value(p_adj))]
  out[]
}

make_hap_test_labels <- function(long, test_dt) {
  if (nrow(test_dt) == 0L) return(data.table::data.table())
  best <- test_dt[!is.na(p_adj)]
  if (nrow(best) == 0L) return(data.table::data.table())
  data.table::setorder(best, trait, p_adj)
  best <- best[, .SD[1L], by = trait]
  y_dt <- long[, .(y = max(value, na.rm = TRUE)), by = trait]
  y_dt[, "y" := y + abs(y) * 0.08 + 1e-8]
  x_dt <- long[, .(x = mean(seq_along(unique(hap_id)))), by = trait]
  merge(best[, .(trait, label)], merge(y_dt, x_dt, by = "trait"), by = "trait")
}

format_p_value <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3g", p)))
}
