demic_env <- new.env(parent = emptyenv())
demic_env$MIN_CONTIGS <- 20
demic_env$MIN_SAMPLES <- 3
demic_env$MAX_ITER <- 3

#' Main function
#'
#' @param X dataframe with coverage matrix
#' (column names: "logCov", "GC", "sample", "contig", "length")
#' @param max_candidate_iter max allowed iterations for estimation of PTR
#' (default: 10)
#' @return dataframe with the estimated PTRs
#' \itemize{
#'   \item estPTR: estimated PTR values
#'   \item coefficient: coefficient of linear regression
#'   \item pValue: p-value of linear regression
#'   \item cor: correlation coefficient
#'   \item correctY: corrected coverage
#' }
#'
#' @examples
#' est_ptrs_001 <- est_ptr(max_bin_001)
#' est_ptrs_001
#'
#' @export
est_ptr <- function(X, max_candidate_iter = 10) {
  verify_input(X)

  contig_est_ptrs <- contigs_pipeline(X)
  sample_est_ptrs <-
    samples_pipeline(X, max_candidate_iter = max_candidate_iter)

  est_ptrs <- combine_ests(contig_est_ptrs, sample_est_ptrs)

  est_ptrs
}
