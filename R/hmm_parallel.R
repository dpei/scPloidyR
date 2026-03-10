#' @title HMM Aggregation Helpers
#' @description Functions for aggregating amplicon-level copy number calls to chromosome level.
#' @name hmm_parallel
NULL

#' @title Calculate Statistical Mode
#' @description Calculates the most frequent value (mode) in a vector.
#' @param x Numeric or integer vector.
#' @param na.rm Logical, whether to remove NA values. Defaults to TRUE.
#' @return The mode value, or NA if vector is empty.
#' @keywords internal
.calculate_mode <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  if (length(x) == 0) {
    return(NA)
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#' @title Aggregate Amplicon-Level Copy Number Calls to Chromosome Level
#' @description Aggregates amplicon-level copy number calls to chromosome-level calls
#'   by calculating the mode of all amplicons within each chromosome.
#' @param hmm_cn_calls Matrix of amplicon-level copy number calls (amplicons x cells).
#' @param read.counts.raw.rowData Data frame containing amplicon metadata with 'chr' column.
#' @return Matrix of chromosome-level copy number calls (chromosomes x cells).
#'
#' @examples
#' \dontrun{
#' # Given amplicon-level CN calls and metadata with chr column
#' cn_by_chr <- aggregate_copy_number(cn_calls_amplicon, amplicon_metadata)
#' }
#' @export
aggregate_copy_number <- function(hmm_cn_calls, read.counts.raw.rowData) {
  if (nrow(hmm_cn_calls) != nrow(read.counts.raw.rowData)) {
    stop("Number of rows in hmm_cn_calls must match number of rows in read.counts.raw.rowData.")
  }

  unique_chrs <- unique(read.counts.raw.rowData$chr)
  cell_barcodes <- colnames(hmm_cn_calls)
  n_cells <- ncol(hmm_cn_calls)

  hmm_cn_by_chr <- matrix(NA, nrow = length(unique_chrs), ncol = n_cells)
  rownames(hmm_cn_by_chr) <- unique_chrs
  colnames(hmm_cn_by_chr) <- cell_barcodes

  for (chrom in unique_chrs) {
    amplicon_indices <- which(read.counts.raw.rowData$chr == chrom)
    chr_amplicon_calls <- hmm_cn_calls[amplicon_indices, , drop = FALSE]
    hmm_cn_by_chr[chrom, ] <- apply(chr_amplicon_calls, 2, .calculate_mode)
  }

  return(hmm_cn_by_chr)
}
