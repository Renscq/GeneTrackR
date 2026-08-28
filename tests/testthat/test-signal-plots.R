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
