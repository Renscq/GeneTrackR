# Author: Rensc
# Date: 2026-05-27
# Version: 0.2.1
# Function: Package-level documentation for GeneTrackR
# Input: GenePred annotations and genomic signal tracks
# Output: Gene model plots, signal plots, and combined genome tracks

#' GeneTrackR: GenePred and Genomic Signal Track Visualization
#'
#' GeneTrackR provides utilities for reading, validating, slicing, merging,
#' writing, and plotting GenePred annotations and bedGraph, wig, or bigWig
#' genomic signal tracks.
#'
#' @keywords internal
#' @import data.table
#' @import ggplot2
#' @importFrom grid arrow unit
#' @importFrom stats median setNames weighted.mean
#' @importFrom utils head modifyList
#' @importFrom rlang .data
#' @importFrom Rcpp evalCpp
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @useDynLib GeneTrackR, .registration = TRUE
"_PACKAGE"

utils::globalVariables(c(
  ".", ":=", "..chrom", "..ext_cols", "..gene_id", "..standard_cols",
  "..strand", "..transcript_id", ".EACHI", ".I", ".N", ".SD", ".data",
  "N", "V1", "V2", "V3", "V4", "bin", "cdsEnd", "cdsEndStat",
  "cdsStart", "cdsStartStat", "cds_end", "cds_end_stat", "cds_start",
  "cds_start_stat", "chrom", "end", "exonCount", "exonFrames", "exon_count",
  "exon_end", "exon_frame", "exon_number", "exon_offset", "exon_start",
  "exon_width", "feature", "gene_end", "gene_id", "gene_start", "gene_type",
  "i.gene_id", "i.sample_id", "i.transcript_id", "mid", "norm_method",
  "old_transcript_id", "reason", "row_id", "sample_id", "score", "source_index",
  "start", "strand",
    "has_strand", "transcript_id", "txEnd", "txStart", "tx_end", "tx_start",
  "length", "length_bp", "value", "width", "x", "xend", "y",
  "line_start", "line_end", "plot_label", "label_x", "label_y",
  "feature_id", "length_plot", "plot_group", "plot_value", "scale_factor", "sample_group",
  "bed_start", "bed_end", "attribute", "phase", "pos", "variant_id",
  "variant_type", "ref", "alt", "qual", "filter", "info", "fill_group",
  "track_source", "label", "ymin", "ymax", "gene_name", "Parent", "Name", "ID",
  "level", "parent_id", "source", "gene_biotype", "transcript_biotype",
  "transcript_type", "fill_group", "rank_pos", "exon_row_key", "strand_key", "phase_int", "cds_start", "cds_end", "tabix_empty_fallback", "tabix_backend", "use_tabix", "has_tabix"
))
