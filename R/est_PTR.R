demic_env <- new.env(parent = emptyenv())
demic_env$MIN_CONTIGS <- 20
demic_env$MIN_SAMPLES <- 3
demic_env$MAX_ITER <- 3

#' Estimate PTRs from both contigs and samples
#'
#' @param X dataframe with coverage matrix
#' (column names: "logCov", "GC", "sample", "contig", "length")
#' @param max_candidate_iter max allowed iterations for estimation of PTR
#' (default: 10)
#' @return named list with results from both methods
#' contigs_ptr dataframe with the estimated PTRs
#' \itemize{
#'   \item estPTR: estimated PTR values
#'   \item coefficient: coefficient of linear regression
#'   \item pValue: p-value of linear regression
#'   \item cor: correlation coefficient
#'   \item correctY: corrected coverage
#' }
#' samples_ptr dataframe with the estimated PTRs
#' \itemize{
#'   \item estPTR: estimated PTR values
#'   \item coefficient: coefficient of linear regression
#'   \item pValue: p-value of linear regression
#'   \item cor: correlation coefficient
#'   \item correctY: corrected coverage
#' }
#'
#' @examples
#' est_ptrs_001 <- est_ptr(ContigCluster1)
#' est_ptrs_001
#'
#' @export
est_ptr <- function(X, max_candidate_iter = 10) {
  verify_input(X)

  contig_est_ptrs <- est_ptr_from_contigs(X)
  sample_est_ptrs <-
    est_ptr_from_samples(X, max_candidate_iter = max_candidate_iter)

  est_ptrs <- combine_ests(contig_est_ptrs, sample_est_ptrs)

  est_ptrs
}
