# Author: Rensc
# Date: 2026-08-30
# Version: dev005
# Function: Package-level documentation for GeneTrackR
# Input: Genomic annotations, signal tracks, variants, and phenotype data
# Output: Integrated genomic tracks, haplotype analyses, and variant prioritization

#' GeneTrackR: Genomic Track Visualization and Haplotype Analysis
#'
#' GeneTrackR provides a shared genomic-coordinate framework for annotation, signal,
#' variant, haplotype, phenotype, linkage-disequilibrium, and variant-effect
#' workflows. Annotation objects interoperate with `GenomicRanges::GRanges`, and
#' bedGraph, wig, and bigWig signal reading, querying, and writing use the built-in
#' pure-R signal I/O layer.
#'
#' @keywords internal
#' @import data.table
#' @importFrom utils head modifyList
#' @importFrom rlang .data
"_PACKAGE"

utils::globalVariables(c(
  ".", ":=", "..ext_cols", "..standard_cols", ".I", ".N", ".SD", ".data", "N", "V1", "V2",
  "V3", "V4", "bin", "cdsEnd", "cdsEndStat", "cdsStart", "cdsStartStat", "cds_end",
  "cds_end_stat", "cds_start", "cds_start_stat", "chrom", "end", "exonCount", "exonFrames",
  "exon_count", "exon_end", "exon_frame", "exon_number", "exon_offset", "exon_start",
  "exon_width", "feature", "gene_end", "gene_id", "gene_start", "gene_type", "i.sample_id",
  "mid", "norm_method", "reason", "row_id", "sample_id", "score", "source_index", "start",
  "strand", "has_strand", "transcript_id", "txEnd", "txStart", "tx_end", "tx_start", "length",
  "length_bp", "value", "width", "x", "xend", "y", "line_start", "line_end", "plot_label",
  "label_x", "label_y", "feature_id", "length_plot", "plot_group", "plot_value",
  "scale_factor", "sample_group", "bed_start", "bed_end", "attribute", "phase", "pos",
  "variant_id", "variant_type", "ref", "alt", "qual", "filter", "info", "fill_group",
  "track_source", "label", "ymin", "ymax", "gene_name", "Parent", "Name", "ID",
  "variant_index", "variant_label", "pair_id", "index_i", "index_j", "variant_i", "variant_j",
  "distance_bp", "n_samples", "D", "Dprime", "Dprime_signed", "p_i", "p_j", "center_x",
  "center_y", "region_x", "genotype_raw", "dosage", "xmin", "xmax", "level", "parent_id",
  "source", "gene_biotype", "transcript_biotype", "transcript_type", "rank_pos",
  "exon_row_key", "strand_key", "phase_int", "tabix_empty_fallback", "tabix_backend",
  "use_tabix", "has_tabix", "tid", "start0", "valid_count", "min_value", "max_value",
  "sum_data", "sum_squared", "size", "chrom_order__", "chrom_size__", "previous_end__"
))
