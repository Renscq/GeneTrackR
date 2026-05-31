# Author: Rensc
# Date: 2026-05-31
# Version: 0.1.0
# Function: Print method for GeneTrackR phenotype plot result
# Input: GeneTrackRPhenoPlot object
# Output: Printed ggplot figure

#' Print GeneTrackR phenotype plot result
#'
#' @param x A GeneTrackRPhenoPlot object returned by plot_hap_pheno() or plot_variant_pheno().
#' @param ... Additional arguments.
#' @return Invisibly returns x.
#' @export
print.GeneTrackRPhenoPlot <- function(x, ...) {
  if (!is.null(x$figure)) {
    print(x$figure)
  }
  invisible(x)
}
