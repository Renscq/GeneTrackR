# Author: Rensc
# Date: 2026-08-31
# Version: dev001
# Function: Package documentation for deterministic GeneTrackR demo files
# Input: None
# Output: R documentation topic

#' GeneTrackR deterministic demo files
#'
#' @description
#' GeneTrackR ships a small deterministic genomic dataset under
#' `inst/extdata`. The files support executable examples, tests, and
#' vignettes without downloading external data.
#'
#' @section Files:
#' The installed `extdata` directory contains:
#' * `gtr_demo.genePredExt`, `gtr_demo.gtf`, and `gtr_demo.gff3` for the same
#'   deterministic gene annotation;
#' * `gtr_demo_features.bed` for genomic interval features;
#' * `gtr_demo_variants.vcf` for designed variant and genotype patterns;
#' * `gtr_demo_pheno.tsv` for deterministic phenotypes;
#' * `gtr_demo.chrom.sizes` for chromosome lengths; and
#' * plus/minus RNA-seq and Ribo-seq `gtr_demo_*bedgraph` signal tracks.
#'
#' @section Provenance:
#' The dataset is synthetic and deterministic. Canonical source tables are
#' stored under `inst/scripts/demo_model` in the source package.
#' `inst/scripts/generate_demo_data.R` generates the installed `gtr_demo_*`
#' files and `inst/scripts/validate_demo_data.R` validates their designed
#' dimensions and relationships. No external biological dataset is bundled.
#'
#' @section Coordinate conventions:
#' Canonical source tables use 1-based closed coordinates. Format-specific
#' output follows the conventions of GenePredExt, BED, bedGraph, GTF, GFF3,
#' and VCF, and GeneTrackR readers convert those inputs to their documented
#' internal coordinate representation.
#'
#' @examples
#' system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR")
#' system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR")
#' system.file("extdata", "gtr_demo_pheno.tsv", package = "GeneTrackR")
#'
#' @name GeneTrackR-demo-data
NULL
