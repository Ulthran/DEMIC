#' Main function
#'
#' @param X dataframe with coverage matrix (column names: "logCov", "GC", "sample", "contig", "length")
#' @param max_candidate_iter max allowed iterations for estimation of PTR (default: 10)
#' @returns dataframe with the estimated PTRs (column names: "estPTR", "coefficient", "pValue", "cor", "correctY")
#'
#' @export
est_ptr <- function(X, max_candidate_iter = 10) {
  est_ptrs <- suppressWarnings({
    contigs_pipeline(X)
  })
  if (is.null(est_ptrs)) {
    Y <- X
    est_ptrs <- samples_pipeline(Y)
  }

  est_ptrs
}
