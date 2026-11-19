#!/usr/bin/env Rscript

# Author: Rensc
# Date: 2026-11-18
# Version: dev001
# Function: Render the compact README and generate its bundled-data figures.
# Input: README.qmd and the current GeneTrackR source checkout.
# Output: README.md and three 500-dpi PNGs under man/figures/.

readme_figure_names <- c(
  "readme-tracks.png", "readme-haplotypes.png", "readme-phenotypes.png"
)

read_readme_source <- function(root) {
  lines <- readLines(file.path(root, "README.qmd"), warn = FALSE, encoding = "UTF-8")
  if (length(lines) == 0L || trimws(lines[1L]) != "---") {
    stop("README.qmd must start with YAML front matter.", call. = FALSE)
  }
  end <- which(trimws(lines[-1L]) == "---")[1L] + 1L
  if (is.na(end) || end >= length(lines)) {
    stop("README.qmd has invalid YAML front matter.", call. = FALSE)
  }
  lines <- lines[seq.int(end + 1L, length(lines))]
  while (length(lines) > 0L && !nzchar(trimws(lines[1L]))) {
    lines <- lines[-1L]
  }
  lines
}

check_readme <- function(root, expected) {
  path <- file.path(root, "README.md")
  if (!file.exists(path) ||
      !identical(readLines(path, warn = FALSE, encoding = "UTF-8"), expected)) {
    stop("README.md is out of sync. Run: Rscript tools/render_readme.R", call. = FALSE)
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("Install the developer dependency 'png' before checking README figures.", call. = FALSE)
  }
  png_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  png_end <- as.raw(c(0, 0, 0, 0, 73, 69, 78, 68, 174, 66, 96, 130))
  for (name in readme_figure_names) {
    image <- file.path(root, "man", "figures", name)
    if (!file.exists(image)) {
      stop(sprintf("Missing README figure: %s", image), call. = FALSE)
    }
    bytes <- readBin(image, what = "raw", n = file.info(image)$size)
    if (length(bytes) < 45L || !identical(bytes[seq_len(8L)], png_signature) ||
        !identical(tail(bytes, 12L), png_end)) {
      stop(sprintf("Invalid README PNG: %s", image), call. = FALSE)
    }
    tryCatch(
      png::readPNG(image, native = TRUE),
      error = function(error) {
        stop(sprintf("Cannot decode README PNG: %s", image), call. = FALSE)
      }
    )
  }
  message("README.md is synchronized and all three PNG assets are present.")
  invisible(TRUE)
}

check_demo <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

readme_text_theme <- function() {
  plain <- ggplot2::element_text(family = "Arial", colour = "black", face = "plain")
  ggplot2::theme(
    text = plain,
    plot.title = plain,
    plot.subtitle = plain,
    plot.caption = plain,
    plot.tag = plain,
    legend.text = plain,
    legend.title = plain,
    strip.text = plain,
    panel.grid = ggplot2::element_blank()
  )
}

generate_readme_figures <- function(root) {
  for (package in c("pkgload", "ragg", "systemfonts")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(sprintf("Install the developer dependency '%s' before generating figures.", package),
           call. = FALSE)
    }
  }
  check_demo(any(systemfonts::system_fonts()$family == "Arial"),
             "Arial must be installed before generating the README figures.")

  # Load this checkout; do not use an unrelated globally installed version.
  pkgload::load_all(root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
  extdata <- function(name) file.path(root, "inst", "extdata", name)
  gp <- GeneTrackR::read_genepred(
    extdata("gtr_demo.genePredExt"), format = "genePredExt",
    verbose = FALSE, progress = FALSE
  )
  signal <- GeneTrackR::read_bwg(
    extdata(c("gtr_demo_rnaseq_plus.bedgraph", "gtr_demo_riboseq_plus.bedgraph")),
    sample_names = c("RNA_seq_plus", "Ribo_seq_plus"), strand = c("+", "+"),
    mode = "memory", verbose = FALSE
  )
  vcf <- GeneTrackR::read_vcf(
    extdata("gtr_demo_variants.vcf"), mode = "memory", keep_genotype = TRUE,
    verbose = FALSE, progress = FALSE
  )
  pheno <- GeneTrackR::read_pheno(
    extdata("gtr_demo_pheno.tsv"), verbose = FALSE, progress = FALSE
  )
  hap <- GeneTrackR::hap_gene_variant(
    vcf, annotation = gp, gene_id = "GeneA", genotype_mode = "string"
  )

  vcf_ids <- as.character(vcf$meta$sample_names)
  hap_ids <- as.character(hap$sample_haplotypes$sample_id)
  phenotype_ids <- as.character(pheno$sample_id)
  check_demo(length(vcf_ids) == 36L && !anyNA(vcf_ids) && !anyDuplicated(vcf_ids),
             "The demo VCF must contain 36 unique nonmissing sample IDs.")
  check_demo(length(phenotype_ids) == 36L && !anyNA(phenotype_ids) &&
               !anyDuplicated(phenotype_ids),
             "The demo phenotype table must contain 36 unique nonmissing sample IDs.")
  check_demo(!anyDuplicated(hap_ids) && setequal(vcf_ids, hap_ids) &&
               setequal(hap_ids, phenotype_ids),
             "VCF, haplotype and phenotype sample IDs must match without record loss.")
  check_demo(nrow(hap$haplotypes) == 4L && all(hap$haplotypes$sample_n == 9L) &&
               nrow(hap$variants) == 11L && !anyDuplicated(hap$variants$variant_id),
             "The GeneA demo must retain 11 variants and four groups of nine samples.")
  check_demo(is.numeric(pheno$protein_content) && all(is.finite(pheno$protein_content)),
             "The demo protein_content trait must contain finite numeric values.")
  check_demo(all(hap$variants$pos >= hap$region$start &
                   hap$variants$pos <= hap$region$end) &&
               !is.unsorted(hap$variants$pos),
             "GeneA variants must be ordered within the 1-based closed gene interval.")
  check_demo(all(signal$data$start >= 1L & signal$data$end >= signal$data$start),
             "Signal intervals must use valid 1-based closed coordinates.")
  message(sprintf(
    "Sample alignment: VCF=%d, phenotype=%d, haplotype=%d, intersection=%d; losses=0/0/0.",
    length(vcf_ids), length(phenotype_ids), length(hap_ids),
    length(intersect(vcf_ids, phenotype_ids))
  ))
  message("GeneA: 11 ordered variants; four haplotypes; n=9 per group; no missing trait values.")

  old_text <- ggplot2::GeomText$default_aes
  on.exit(ggplot2::update_geom_defaults("text", old_text), add = TRUE)
  ggplot2::update_geom_defaults(
    "text", list(family = "Arial", colour = "black", fontface = "plain")
  )
  set.seed(1014)
  tracks <- GeneTrackR::plot_tracks(
    annotation = gp, signal = signal, gene_id = "GeneA", strand = "+",
    signal_y_ticks = "pretty",
    signal_colors = c(RNA_seq_plus = "#4477AA", Ribo_seq_plus = "#228833"),
    layout = "gene_top", plot_theme = "classic", text_size = 11
  ) & readme_text_theme()
  haplotypes <- GeneTrackR::plot_hap_variant(
    hap, annotation = gp, min_hap_samples = 3, variant_label = "pos",
    gene_track_legend_position = "top", text_size = 11,
    genotype_text_size = 3.2, table_x_angle = 90
  ) & (readme_text_theme() +
         ggplot2::theme(plot.margin = ggplot2::margin(6, 30, 6, 6)))
  association <- GeneTrackR::plot_hap_pheno(
    hap, phenotype = pheno, traits = "protein_content", min_hap_samples = 3,
    test_method = "t.test", p_adjust = "BH", p_value_type = "adjusted",
    show_points = TRUE, x_text_angle = 0, text_size = 11,
    fill_colors = c(Hap1 = "#4477AA", Hap2 = "#66CCEE", Hap3 = "#228833", Hap4 = "#CCBB44")
  )
  check_demo(nrow(association$plot_data) == 36L && nrow(association$pvalue) == 6L,
             "The phenotype plot must retain all 36 samples and six pairwise tests.")
  phenotype_plot <- association$figure +
    ggplot2::theme_classic(base_family = "Arial", base_size = 11) +
    readme_text_theme() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(colour = "black", face = "plain"),
      axis.title = ggplot2::element_text(colour = "black", face = "plain")
    ) +
    ggplot2::labs(x = "GeneA haplotype", y = "Protein content (arbitrary units)", title = NULL)

  figure_dir <- file.path(root, "man", "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  plots <- list(tracks, haplotypes, phenotype_plot)
  widths <- c(11, 11, 7)
  heights <- c(5.5, 6.5, 5)
  for (i in seq_along(plots)) {
    ggplot2::ggsave(
      filename = file.path(figure_dir, readme_figure_names[i]), plot = plots[[i]],
      device = ragg::agg_png, width = widths[i], height = heights[i], units = "in",
      dpi = 500, bg = "white"
    )
    message(sprintf("Generated man/figures/%s at 500 dpi.", readme_figure_names[i]))
  }
  invisible(association)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (identical(args, "--help")) {
    cat("Usage: Rscript tools/render_readme.R [--check | --help]\n",
        "Run from the package root. Default: generate real figures and README.md.\n",
        "--check validates Markdown synchronization and committed PNG assets without plotting.\n",
        "README.qmd is used as Markdown source; this script does not invoke Quarto.\n",
        sep = "")
    return(invisible(NULL))
  }
  if (length(args) > 1L || (length(args) == 1L && args != "--check")) {
    stop("Usage: Rscript tools/render_readme.R [--check | --help]", call. = FALSE)
  }
  root <- normalizePath(".", winslash = "/", mustWork = TRUE)
  expected <- read_readme_source(root)
  if (!identical(args, "--check")) {
    generate_readme_figures(root)
    writeLines(expected, file.path(root, "README.md"), useBytes = TRUE)
    message("Rendered README.md from README.qmd.")
  }
  check_readme(root, expected)
}

if (sys.nframe() == 0L) {
  tryCatch(main(), error = function(error) {
    message(conditionMessage(error))
    quit(status = 1L)
  })
}
