test_that("read_bwg reads strand-specific demo bedGraph tracks", {
  files <- gtr_extdata(c(
    "gtr_demo_rnaseq_plus.bedgraph",
    "gtr_demo_rnaseq_minus.bedgraph",
    "gtr_demo_riboseq_plus.bedgraph",
    "gtr_demo_riboseq_minus.bedgraph"
  ))
  bg <- read_bwg(
    files,
    format = "bedgraph",
    sample_names = c("RNA_seq_plus", "RNA_seq_minus", "Ribo_seq_plus", "Ribo_seq_minus"),
    strand = c("+", "-", "+", "-"),
    mode = "lazy",
    verbose = FALSE
  )

  expect_s3_class(bg, "BwgTrack")
  expect_equal(nrow(bg$samples), 4L)
  expect_identical(bg$samples$strand, c("+", "-", "+", "-"))

  dt_plus <- retrieve_bwg(
    bg,
    chrom = "chr1",
    start = 12340001,
    end = 12341000,
    strand = "+"
  )
  expect_gt(nrow(dt_plus), 0L)
  expect_true(all(dt_plus$strand == "+"))

  dt_minus <- retrieve_bwg(
    bg,
    chrom = "chr1",
    start = 12356001,
    end = 12357500,
    strand = "-"
  )
  expect_gt(nrow(dt_minus), 0L)
  expect_true(all(dt_minus$strand == "-"))
})

test_that("gene model, signal, and combined track plots use unified demo data", {
  gp <- read_genepred(
    gtr_extdata("gtr_demo.genePredExt"),
    format = "genePredExt",
    verbose = FALSE
  )
  signal_all <- read_bwg(
    gtr_extdata(c(
      "gtr_demo_rnaseq_plus.bedgraph",
      "gtr_demo_rnaseq_minus.bedgraph",
      "gtr_demo_riboseq_plus.bedgraph",
      "gtr_demo_riboseq_minus.bedgraph"
    )),
    format = "bedgraph",
    sample_names = c("RNA_seq_plus", "RNA_seq_minus", "Ribo_seq_plus", "Ribo_seq_minus"),
    strand = c("+", "-", "+", "-"),
    mode = "lazy",
    verbose = FALSE
  )
  bed <- read_bed(gtr_extdata("gtr_demo_features.bed"), verbose = FALSE)
  vcf <- read_vcf(gtr_extdata("gtr_demo_variants.vcf"), mode = "memory", verbose = FALSE)

  expect_s3_class(plot_gene(gp, gene_id = "GeneA"), "ggplot")
  expect_s3_class(plot_transcript(gp, transcript_id = "TxA1"), "ggplot")
  expect_s3_class(
    plot_region(gp, chrom = "chr1", start = 12339001, end = 12374500),
    "ggplot"
  )

  expect_true(inherits(plot_signal_gene(signal_all, gp, gene_id = "GeneA"), "patchwork"))
  expect_true(inherits(plot_signal_transcript(signal_all, gp, transcript_id = "TxA1"), "patchwork"))
  expect_true(inherits(
    plot_signal_region(
      signal_all,
      chrom = "chr1",
      start = 12339001,
      end = 12374500,
      strand = "both",
      annotation = gp
    ),
    "patchwork"
  ))

  p <- plot_tracks(
    annotation = gp,
    signal = signal_all,
    features = bed,
    variants = vcf,
    chrom = "chr1",
    start = 12339001,
    end = 12374500,
    strand = "both",
    gene_palette = "Set2"
  )
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("signal plotting supports assay groups and summaries", {
  gp <- read_genepred(
    gtr_extdata("gtr_demo.genePredExt"),
    format = "genePredExt",
    verbose = FALSE
  )
  signal_all <- read_bwg(
    gtr_extdata(c(
      "gtr_demo_rnaseq_plus.bedgraph",
      "gtr_demo_rnaseq_minus.bedgraph",
      "gtr_demo_riboseq_plus.bedgraph",
      "gtr_demo_riboseq_minus.bedgraph"
    )),
    format = "bedgraph",
    sample_names = c("RNA_seq_plus", "RNA_seq_minus", "Ribo_seq_plus", "Ribo_seq_minus"),
    strand = c("+", "-", "+", "-"),
    mode = "lazy",
    verbose = FALSE
  )
  groups <- c(
    RNA_seq_plus = "RNA-seq",
    RNA_seq_minus = "RNA-seq",
    Ribo_seq_plus = "Ribo-seq",
    Ribo_seq_minus = "Ribo-seq"
  )

  p <- plot_signal_region(
    signal_all,
    annotation = gp,
    chrom = "chr1",
    start = 12339001,
    end = 12374500,
    strand = "both",
    sample_groups = groups,
    signal_color_by = "group",
    signal_summary = "mean",
    bin_size = 100
  )
  expect_true(inherits(p, "patchwork"))
})

test_that("track themes and signal style parameters are applied consistently", {
  dt <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = c(1L, 3L),
    end = c(2L, 4L),
    value = c(4, 25),
    strand = "*"
  )

  p <- plot_signal_core(
    dt,
    plot_type = "bar",
    signal_transform = "sqrt",
    signal_y_scale = "fixed",
    signal_y_limits = c(0, 4),
    signal_alpha = 0.60,
    signal_bar_width = 0.50,
    plot_theme = "classic",
    show_panel_border = FALSE
  )
  built <- ggplot2::ggplot_build(p)
  bar_layer <- built$data[[1L]]

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(bar_layer), 2L)
  expect_equal(bar_layer$xmax - bar_layer$xmin, rep(1, 2))
  expect_equal(unique(bar_layer$alpha), 0.60)
  expect_s3_class(p$theme$panel.border, "element_blank")
  expect_equal(p$scales$get_scales("y")$limits, c(0, 4))

  expect_warning(
    plot_signal_core(
      dt,
      signal_y_scale = "free",
      signal_y_limits = c(0, 10)
    ),
    "changed to 'fixed'"
  )
})

test_that("new plotting parameter validation is strict", {
  expect_error(normalize_signal_alpha(1.1), "between 0 and 1")
  expect_error(normalize_signal_bar_width(0), "greater than 0")
  expect_error(normalize_signal_y_limits(c(2, 1)), "increasing")
  expect_error(normalize_show_panel_border(1), "NULL, TRUE, or FALSE")
  expect_error(normalize_plot_theme("dark"), "arg")
})

test_that("gene model feature colors are stable across observed feature subsets", {
  utr_only <- make_gene_model_fill_colors(features = "UTR")
  cds_utr <- make_gene_model_fill_colors(features = c("CDS", "UTR"))
  exon_only <- make_gene_model_fill_colors(features = "exon")

  expect_identical(gene_model_feature_levels(), c("UTR", "CDS", "exon"))
  expect_true(all(c("UTR", "CDS", "exon") %in% names(utr_only)))
  expect_identical(names(utr_only)[seq_along(gene_model_feature_levels())], gene_model_feature_levels())
  expect_identical(unname(utr_only["UTR"]), unname(cds_utr["UTR"]))
  expect_identical(unname(exon_only["exon"]), unname(cds_utr["exon"]))

  alias_cols <- make_gene_model_fill_colors(features = c("five_prime_UTR", "three_prime_UTR"))
  expect_true("UTR" %in% names(alias_cols))
})


test_that("variant marker colors are stable and ordered", {
  cols <- make_variant_marker_fill_colors(variant_types = c("DEL", "SNP"))
  expect_identical(variant_marker_levels(), c("SNP", "Ind", "..."))
  expect_true(all(c("SNP", "Ind", "...") %in% names(cols)))
  expect_identical(normalize_variant_marker_type(c("SNP", "INS", "unknown")), c("SNP", "Ind", "..."))
})

test_that("discrete signal palettes preserve standard class order", {
  blues3 <- RColorBrewer::brewer.pal(3L, "Blues")
  blues4 <- RColorBrewer::brewer.pal(4L, "Blues")

  expect_identical(
    make_signal_palette(2L, "Blues", signal_palette_direction = 1),
    unname(blues3[1:2])
  )
  expect_identical(
    make_signal_palette(2L, "Blues", signal_palette_direction = -1),
    unname(rev(blues3)[1:2])
  )
  expect_identical(
    make_signal_palette(4L, "Blues", signal_palette_direction = 1),
    unname(blues4)
  )

  continuous <- make_signal_continuous_palette(
    32L,
    "Blues",
    signal_palette_direction = 1
  )
  continuous_rev <- make_signal_continuous_palette(
    32L,
    "Blues",
    signal_palette_direction = -1
  )
  expect_identical(continuous, rev(continuous_rev))

  exact <- c(sampleA = "#112233", sampleB = "#445566")
  expect_identical(
    normalize_signal_colors(
      c("sampleA", "sampleB"),
      signal_palette = "Blues",
      signal_palette_direction = -1,
      signal_colors = exact
    ),
    exact
  )
})



test_that("qualitative signal palettes preserve palette order", {
  set1 <- RColorBrewer::brewer.pal(3L, "Set1")
  expected <- stats::setNames(set1, c("frame0", "frame1", "frame2"))

  expect_identical(
    make_signal_palette(3L, "Set1", signal_palette_direction = 1),
    unname(set1)
  )
  expect_identical(
    make_signal_palette(3L, "Set1", signal_palette_direction = -1),
    rev(unname(set1))
  )
  expect_identical(
    make_frame_colors(frame_palette = "Set1"),
    expected
  )

  ordered_ids <- c("sampleB", "sampleA", "sampleC")
  ordered_cols <- normalize_signal_colors(
    ordered_ids,
    signal_palette = "Set1"
  )
  expect_identical(names(ordered_cols), ordered_ids)
  expect_identical(unname(ordered_cols), unname(set1))
})

test_that("plot_signal_region uses sample order for discrete palette mapping", {
  sample_ids <- c("RNA_seq_plus", "RNA_seq_minus", "Ribo_seq_plus", "Ribo_seq_minus")
  expected <- stats::setNames(
    RColorBrewer::brewer.pal(4L, "Set1"),
    sample_ids
  )
  observed <- normalize_signal_colors(
    sample_ids,
    signal_palette = "Set1",
    signal_palette_direction = 1
  )

  expect_identical(observed, expected)
})


test_that("signal group palette order follows sample-group mapping order", {
  dt <- data.table::data.table(
    sample_id = factor(
      c("sampleB", "sampleA", "sampleC"),
      levels = c("sampleA", "sampleB", "sampleC")
    ),
    chrom = "chr1",
    start = 1:3,
    end = 1:3,
    value = 1:3,
    strand = "*"
  )
  groups <- c(sampleA = "group2", sampleB = "group1", sampleC = "group3")
  grouped <- apply_signal_grouping(dt, sample_groups = groups, signal_summary = "none")

  expect_true(is.factor(grouped$sample_group))
  expect_identical(levels(grouped$sample_group), c("group2", "group1", "group3"))
})

test_that("signal plot wrappers accept explicit signal and gene track heights", {
  gp <- read_genepred(
    gtr_extdata("gtr_demo.genePredExt"),
    format = "genePredExt",
    verbose = FALSE
  )
  rnaseq <- read_bwg(
    gtr_extdata(c(
      "gtr_demo_rnaseq_plus.bedgraph",
      "gtr_demo_rnaseq_minus.bedgraph"
    )),
    format = "bedgraph",
    sample_names = c("RNA_seq_plus", "RNA_seq_minus"),
    strand = c("+", "-"),
    mode = "lazy",
    verbose = FALSE
  )

  expect_true(inherits(
    plot_signal_gene(
      rnaseq,
      gp,
      gene_id = "GeneA",
      signal_palette_direction = -1,
      signal_track_height = 4,
      gene_track_height = 1
    ),
    "patchwork"
  ))
  expect_true(inherits(
    plot_signal_transcript(
      rnaseq,
      gp,
      transcript_id = "TxA1",
      signal_track_height = 2,
      gene_track_height = 1
    ),
    "patchwork"
  ))
  expect_true(inherits(
    plot_signal_region(
      rnaseq,
      chrom = "chr1",
      start = 12339001,
      end = 12374500,
      annotation = gp,
      strand = "both",
      signal_track_height = 5,
      gene_track_height = 1
    ),
    "patchwork"
  ))

  expect_equal(normalize_track_height(4, "signal_track_height", 3), 4)
  expect_warning(
    expect_equal(normalize_track_height(0, "signal_track_height", 3), 3),
    "positive numeric"
  )
})


test_that("feature browser tracks use compact automatic legend groups", {
  bed <- read_bed(gtr_extdata("gtr_demo_features.bed"), verbose = FALSE)
  dt <- retrieve_feature(
    bed,
    chrom = "chr1",
    start = 12339001,
    end = 12374500,
    mode = "trim",
    level = "feature",
    as = "data.table"
  )

  expect_identical(
    extract_feature_group_from_name(c("GeneA_promoter|promoter", "QTL1|QTL", "plain_name")),
    c("promoter", "QTL", NA_character_)
  )

  expect_identical(
    resolve_feature_color_column(dt, color_by = "auto"),
    "feature_group"
  )

  compact <- compact_feature_levels(
    c("promoter", "enhancer", "candidate_region", "QTL", "repeat", "conserved_region"),
    max_legend_levels = 5
  )
  expect_lte(length(levels(compact)), 5L)

  p <- plot_feature_track(
    bed,
    chrom = "chr1",
    start = 12339001,
    end = 12374500,
    color_by = "auto",
    max_legend_levels = 5
  )
  expect_s3_class(p, "ggplot")
})

test_that("exported plotting palette defaults follow package conventions", {
  plot_functions <- grep(
    "^plot_",
    getNamespaceExports("GeneTrackR"),
    value = TRUE
  )

  for (fn_name in plot_functions) {
    fn <- getExportedValue("GeneTrackR", fn_name)
    defaults <- formals(fn)
    palette_args <- grep("palette$", names(defaults), value = TRUE)
    for (arg_name in palette_args) {
      expected_palette <- if (
        identical(fn_name, "plot_ld_block") &&
          identical(arg_name, "color_palette")
      ) {
        "Reds"
      } else {
        "Paired"
      }
      expect_identical(
        defaults[[arg_name]],
        expected_palette,
        info = paste(
          fn_name,
          arg_name,
          "should default to",
          expected_palette
        )
      )
    }
  }
})

test_that("plot_tracks forwards variant palette and explicit colors", {
  variants <- read_vcf(
    gtr_extdata("gtr_demo_variants.vcf"),
    mode = "memory",
    verbose = FALSE,
    progress = FALSE
  )
  variant_colors <- c(
    SNP = "#1F78B4",
    INS = "#33A02C",
    DEL = "#E31A1C"
  )

  plots <- make_variant_track_plots(
    variants,
    chrom = "chr1",
    start = 12339001,
    end = 12374500,
    variant_palette = "Paired",
    variant_colors = variant_colors
  )
  expect_length(plots, 1L)
  expect_s3_class(plots[[1L]], "ggplot")

  built <- ggplot2::ggplot_build(plots[[1L]])
  observed <- unique(unlist(lapply(built$data, function(x) x$colour), use.names = FALSE))
  observed <- observed[!is.na(observed)]
  expect_true(all(observed %in% unname(variant_colors)))
})

signal_plot_strand_fixture <- function() {
  BwgTrack(
    samples = data.frame(
      sample_id = c("plus", "minus"), strand = c("+", "-"), has_strand = TRUE
    ),
    data = data.frame(
      sample_id = c("plus", "minus", "minus", "minus"), chrom = "chr1",
      start = c(25L, 21L, 24L, 30L), end = c(29L, 22L, 29L, 32L),
      value = c(99, 5, 2, 7), strand = c("+", "-", "-", "-")
    )
  )
}

test_that("negative strand signal uses genomic gene and spliced transcript coordinates", {
  track <- signal_plot_strand_fixture()
  annotation <- coverage_gene_model()
  gene <- plot_signal_gene(track, annotation, "minus_gene", show_gene_model = FALSE)
  transcript <- plot_signal_transcript(
    track, annotation, "minus_tx", show_gene_model = FALSE
  )

  expect_identical(levels(gene$data$sample_id), "minus")
  expect_equal(gene$data$start, c(21L, 24L, 30L))
  expect_equal(gene$data$end, c(22L, 29L, 32L))
  expect_equal(gene$data$value, c(5, 2, 7))
  gene_plot <- ggplot2::ggplot_build(gene)
  expect_equal(gene_plot$layout$panel_params[[1L]]$x.range, c(20.45, 32.55))
  gene_bars <- gene_plot$data[[2L]]
  expect_equal(gene_bars$xmin, c(20.5, 23.5, 29.5))
  expect_equal(gene_bars$xmax, c(22.5, 29.5, 32.5))

  expect_identical(levels(transcript$data$sample_id), "minus")
  # The interval from 24 to 29 spans the intron and contributes to both exons.
  expect_equal(transcript$data$start, c(1L, 4L, 6L, 7L))
  expect_equal(transcript$data$end, c(2L, 5L, 6L, 9L))
  expect_equal(transcript$data$value, c(5, 2, 2, 7))
  tx_plot <- ggplot2::ggplot_build(transcript)
  expect_equal(tx_plot$layout$panel_params[[1L]]$x.range, c(0.6, 9.4))
  tx_bars <- tx_plot$data[[2L]]
  expect_equal(tx_bars$xmin, c(0.5, 3.5, 5.5, 6.5))
  expect_equal(tx_bars$xmax, c(2.5, 5.5, 6.5, 9.5))
})

test_that("frame plots count negative strand CDS bases in biological order", {
  track <- signal_plot_strand_fixture()
  annotation <- coverage_gene_model()
  frame_colors <- c(frame0 = "#0072B2", frame1 = "#D55E00", frame2 = "#009E73")
  expect_message(
    transcript <- plot_signal_transcript(
      track, annotation, "minus_tx", plot_type = "frame", bin_size = 3L,
      frame_colors = frame_colors[c(3L, 1L, 2L)],
      signal_y_ticks = "pretty", show_gene_model = FALSE
    ),
    "requires base-resolution signal"
  )
  layers <- ggplot2::ggplot_build(transcript)$data
  utr <- layers[[1L]]
  cds <- layers[[2L]]
  cds <- cds[order(cds$x), ]
  # CDS bases are 30, 29, 25, 24, 23 in biological order.
  expected_colors <- unname(frame_colors[c(1L, 3L, 2L, 1L)])
  expect_equal(utr$x, c(1, 2, 8, 9))
  expect_equal(utr$y, c(5, 5, 7, 7))
  expect_equal(utr$xmin, c(0.5, 1.5, 7.5, 8.5))
  expect_equal(utr$xmax, c(1.5, 2.5, 8.5, 9.5))
  expect_false(any(utr$fill %in% frame_colors))
  expect_equal(cds$x, c(4, 5, 6, 7))
  expect_equal(cds$y, c(2, 2, 2, 7))
  expect_identical(cds$fill, expected_colors)

  genomic <- plot_signal_transcript(
    track, annotation, "minus_tx", coordinate = "genomic", plot_type = "frame",
    frame_colors = frame_colors, signal_y_ticks = "pretty", show_gene_model = FALSE
  )
  genomic_cds <- ggplot2::layer_data(genomic, 2L)
  genomic_cds <- genomic_cds[order(genomic_cds$x), ]
  expect_equal(genomic_cds$x, c(24, 25, 29, 30))
  expect_identical(genomic_cds$fill, expected_colors)

  intron <- slice_bwg(track, "chr1", 26L, 28L, samples = "minus")
  expect_error(
    plot_signal_transcript(intron, annotation, "minus_tx", plot_type = "frame"),
    "No signal records overlapped the exon positions"
  )
})

test_that("named sample groups and signal colors follow names in any order", {
  track <- signal_plot_strand_fixture()
  groups <- c(plus = "treated", minus = "control")
  colors <- c(treated = "#D55E00", control = "#0072B2")
  expected_groups <- c("treated", rep("control", 3L))

  for (index in list(1:2, 2:1)) {
    plot <- plot_signal_region(
      track, "chr1", 21L, 32L, strand = "both", sample_groups = groups[index],
      signal_color_by = "group", signal_colors = colors[rev(index)],
      signal_palette_direction = -1, signal_y_scale = "fixed"
    )
    rows <- plot$data[order(plot$data$sample_id, plot$data$start)]
    expect_identical(as.character(rows$sample_group), expected_groups)
    bars <- ggplot2::layer_data(plot)
    bars <- bars[order(bars$PANEL, bars$xmin), ]
    expect_identical(bars$fill, unname(colors[expected_groups]))
  }
})

test_that("signal plots aggregate normalized replicates within bins before transformation", {
  track <- BwgTrack(
    samples = data.frame(sample_id = c("reference", "rep_b", "rep_a", "rep_c")),
    data = data.frame(
      sample_id = rep(c("rep_a", "rep_b", "rep_c", "reference"), each = 3L),
      chrom = "chr1", strand = "*", start = rep(c(2L, 4L, 8L), 4L),
      end = rep(c(3L, 5L, 9L), 4L),
      value = c(-6, -3, 3, 1, 3, 2, 3, 5, 4, 1, 3, 0)
    )
  )
  raw <- as.data.frame(track$data)
  normalized <- norm_bwg(track, method = "scale", scale_factor = c(rep_a = 2, reference = 0.5))
  groups <- data.frame(sample_id = c("rep_c", "rep_a", "rep_b"), group = "replicates")
  group_values <- list(mean = c(-1, 4), median = c(2, 4), sum = c(-3, 12))
  colors <- c("#0072B2", "#D55E00")

  for (summary in names(group_values)) {
    plot <- plot_signal_region(
      normalized, "chr1", 1L, 10L, bin_size = 4L, sample_groups = groups,
      signal_summary = summary, signal_color_by = "group", signal_colors = colors,
      signal_transform = "log2", signal_y_scale = "fixed"
    )
    rows <- plot$data[order(plot$data$sample_id, plot$data$start)]
    expected <- c(1, 0, group_values[[summary]])
    transformed <- sign(expected) * log2(abs(expected) + 1)
    expect_identical(as.character(rows$sample_id), rep(c("reference", "replicates"), each = 2L))
    expect_equal(rows$start, rep(c(2L, 8L), 2L))
    expect_equal(rows$end, rep(c(5L, 9L), 2L))
    expect_equal(rows$value, expected)
    expect_equal(rows$plot_value, transformed)

    bars <- ggplot2::layer_data(plot)
    bars <- bars[order(bars$PANEL, bars$xmin), ]
    expect_equal(bars$ymin, pmin(transformed, 0))
    expect_equal(bars$ymax, pmax(transformed, 0))
    expect_identical(bars$fill, rep(colors, each = 2L))
  }
  expect_identical(as.data.frame(track$data), raw)

  expect_error(
    plot_signal_region(track, "chr1", 1L, 10L, sample_groups = data.frame(sample_id = "rep_a")),
    "must contain `sample_id` and `group` columns"
  )
  expect_error(
    plot_signal_region(track, "chr1", 1L, 10L, sample_groups = "replicates"),
    "must have the same length as selected samples"
  )
})

test_that("line signal returns to zero across gaps and preserves empty selected samples", {
  track <- BwgTrack(
    samples = data.frame(sample_id = c("empty", "covered", "excluded")),
    data = data.frame(
      sample_id = c("covered", "covered", "empty", "excluded"),
      chrom = c("chr1", "chr1", "chr2", "chr1"), strand = "*",
      start = c(2L, 7L, 5L, 5L), end = c(3L, 8L, 5L, 5L),
      value = c(-9, 99, 7, 10)
    )
  )
  selected <- c("empty", "covered")
  plot <- plot_signal_region(
    track, "chr1", 1L, 10L, samples = selected, plot_type = "line",
    signal_transform = "log10", signal_y_ticks = "pretty",
    sample_groups = c("shared", "shared"), signal_color_by = "group"
  )
  built <- ggplot2::ggplot_build(plot)
  panels <- built$layout$layout
  expect_identical(as.character(panels$sample_id), selected)
  line <- built$data[[1L]]
  ids <- panels$sample_id[match(line$PANEL, panels$PANEL)]
  covered <- line[ids == "covered", ]
  expect_equal(covered$x, c(2.5, 3.5, 6.5, 7.5))
  expect_equal(covered$y, c(-1, 0, 0, 2))
  expect_length(unique(covered$group), 1L)
  expect_equal(line$y[ids == "empty"], 0)

  absent <- plot_signal_region(track, "chr3", 1L, 10L, samples = selected)$data
  expect_identical(as.character(absent$sample_id), selected)
  expect_equal(absent$start, c(1L, 1L))
  expect_equal(absent$end, c(10L, 10L))
  expect_equal(absent$value, c(0, 0))
  expect_error(
    plot_signal_region(track, "chr1", 1L, 10L, samples = "unknown"),
    "No signal records were found"
  )
})

test_that("heatmap signal preserves sample order and explicit or adaptive bin widths", {
  positions <- rep(c(1L, 2L, 5L, 6L, 9L), 2L)
  track <- BwgTrack(
    samples = data.frame(sample_id = c("B", "A")),
    data = data.frame(
      sample_id = rep(c("A", "B"), each = 5L), chrom = "chr1", strand = "*",
      start = positions, end = positions, value = c(3, 15, 1, 3, 7, -3, -15, 1, 0, 15)
    )
  )
  colors <- c("#2166AC", "#B2182B")
  plot <- plot_signal_region(
    track, "chr1", 1L, 9L, plot_type = "heatmap", heatmap_bin_size = 4L,
    heatmap_summary = "sum", signal_colors = colors
  )
  rows <- plot$data[order(plot$data$sample_id, plot$data$start)]
  expect_identical(as.character(rows$sample_id), rep(c("B", "A"), each = 3L))
  expect_equal(rows$start, rep(c(1L, 5L, 9L), 2L))
  expect_equal(rows$end, rep(c(4L, 8L, 9L), 2L))
  expect_equal(rows$value, c(-18, 1, 15, 18, 4, 7))

  tiles <- ggplot2::layer_data(plot)
  tiles <- tiles[order(tiles$y, tiles$x), ]
  expect_equal(as.numeric(tiles$y), rep(1:2, each = 3L))
  expect_equal(tiles$xmax - tiles$xmin, rep(c(4, 4, 1), 2L))
  expect_identical(tiles$fill[c(1L, 4L)], colors)

  adaptive <- plot_signal_region(
    track, "chr1", 1L, 9L, plot_type = "heatmap", heatmap_max_bins = 3L,
    heatmap_summary = "max", signal_transform = "log2"
  )$data
  adaptive <- adaptive[order(adaptive$sample_id, adaptive$start)]
  expect_equal(adaptive$start, rep(c(1L, 4L, 7L), 2L))
  expect_equal(adaptive$end, rep(c(3L, 6L, 9L), 2L))
  expect_equal(adaptive$value, c(-3, 1, 15, 15, 3, 7))
  expect_equal(adaptive$plot_value, c(-2, 1, 4, 4, 2, 3))
})
