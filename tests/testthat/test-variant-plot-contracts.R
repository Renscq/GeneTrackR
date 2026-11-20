variant_plot_fixture <- function() {
    VariantTrack(data.frame(
        chrom = c("chr1", "chr1", "chr2"), pos = c(5L, 10L, 3L),
        variant_id = c("snp1", "ins1", "snp2"), ref = c("A", "A", "C"),
        alt = c("G", "AT", "T"), filter = c("PASS", NA, "q10")
    ))
}

test_that("variant plots preserve regional selection and explicit label meanings", {
    variant <- variant_plot_fixture()
    before <- serialize(variant, NULL)
    for (label in c("variant_id", "variant_type", "ref", "alt")) {
        plot <- plot_variant(variant, chrom = "chr1", start = 5L, end = 10L,
            color_by = "filter", label_by = label, variant_shape = "point",
            variant_colors = c("#0072B2", "#D55E00"), plot_theme = "classic")
        expect_identical(plot$data$variant_id, c("snp1", "ins1"))
        expect_identical(plot$data$label, as.character(plot$data[[label]]))
        expect_identical(plot$data$color_group, c("PASS", "NA"))
        built <- ggplot2::ggplot_build(plot)
        expect_equal(built$data[[1L]]$x, c(5, 10))
        expect_identical(built$data[[1L]]$colour, c("#0072B2", "#D55E00"))
        expect_identical(built$data[[2L]]$label, plot$data$label)
    }
    expect_identical(serialize(variant, NULL), before)
})

test_that("variant shape and palette choices retain positions and alias output", {
    variant <- variant_plot_fixture()
    for (shape in c("lollipop", "point", "rug")) {
        plot <- plot_variant(variant, chrom = "chr1", variant_type = "SNP",
            color_by = "none", variant_shape = shape)
        built <- ggplot2::ggplot_build(plot)
        expect_equal(built$data[[1L]]$x, 5)
        expect_identical(plot$data$color_group, "variant")
        expect_equal(length(plot$layers), if (shape == "lollipop") 2L else 1L)
    }
    original <- plot_variant(variant, "chr1", 1L, 20L, label_by = "variant_id")
    alias <- plot_variant_track(variant, "chr1", 1L, 20L, label_by = "variant_id")
    expect_equal(alias$data, original$data)
    expect_equal(ggplot2::ggplot_build(alias)$data,
        ggplot2::ggplot_build(original)$data)
    fallback <- plot_variant(variant, variant_palette = "not-a-brewer-palette")
    expect_equal(nrow(ggplot2::ggplot_build(fallback)$data[[1L]]), 3L)
})

test_that("empty variant plots retain requested axes and reject invalid regions", {
    variant <- variant_plot_fixture()
    plot <- plot_variant(variant, chrom = "chr9", start = 20L, end = 30L,
        plot_theme = "classic")
    expect_equal(plot$coordinates$limits$x, c(20L, 30L))
    expect_identical(plot$labels$x, "Chromosome chr9 position (bp)")
    expect_equal(nrow(ggplot2::ggplot_build(plot)$data[[1L]]), 2L)
    empty <- VariantTrack(variant$data[0L])
    inferred <- plot_variant(empty)
    expect_equal(inferred$coordinates$limits$x, c(0L, 1L))
    expect_identical(inferred$labels$x, "Genomic position (bp)")
    expect_error(plot_variant(variant, "chr1", 20L, 10L), "start|region")
    expect_error(plot_variant(list()), "VariantTrack")
})
