#' Main function
#'
#' @param X dataframe with coverage matrix (column names: "logCov", "GC", "sample", "contig", "length")
#' @param max_candidate_iter max allowed iterations for estimation of PTR (default: 10)
#' @returns dataframe with the estimated PTRs (column names: "estPTR", "coefficient", "pValue", "cor", "correctY")
#'
#' @export
est_PTR <- function(X, max_candidate_iter = 10) {
  est_PTRs <- suppressWarnings({
    contigs_pipeline(X)
  })
  if (is.null(est_PTRs)) {
    Y <- X
    est_PTRs <- samples_pipeline(Y)
  }

  est_PTRs
}
