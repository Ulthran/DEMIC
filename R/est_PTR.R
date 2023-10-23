#' Main function
#'
#' @param X dataframe with coverage matrix
#' (column names: "logCov", "GC", "sample", "contig", "length")
#' @param max_candidate_iter max allowed iterations for estimation of PTR
#' (default: 10)
#' @returns dataframe with the estimated PTRs
#' \itemize{
#'   \item estPTR: estimated PTR values
#'   \item coefficient: coefficient of linear regression
#'   \item pValue: p-value of linear regression
#'   \item cor: correlation coefficient
#'   \item correctY: corrected coverage
#' }
#'
#' @export
est_ptr <- function(X, max_candidate_iter = 10) {
  est_ptrs <- suppressWarnings({
    contigs_pipeline(X)
  })
  if (is.null(est_ptrs)) {
    Y <- X
    browser()
    est_ptrs <- samples_pipeline(Y)
  }

  est_ptrs
}
