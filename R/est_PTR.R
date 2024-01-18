demic_env <- new.env(parent = emptyenv())
demic_env$MIN_CONTIGS <- 20
demic_env$MIN_SAMPLES <- 3
demic_env$MAX_ITER <- 3

#' Estimate PTRs from both contigs and samples
#'
#' @param X dataframe with coverage matrix
#' (column names: "log_cov", "GC_content", "sample", "contig", "length")
#' @return named list with results from both methods
#' contigs_ptr dataframe with the estimated PTRs
#' \itemize{
#'   \item est_ptr: estimated PTR values
#'   \item coefficient: coefficient of linear regression
#'   \item pValue: p-value of linear regression
#'   \item cor: correlation coefficient
#'   \item correctY: corrected coverage
#' }
#' samples_ptr dataframe with the estimated PTRs
#' \itemize{
#'   \item est_ptr: estimated PTR values
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
est_ptr <- function(X) {
  verify_input(X)

  tryCatch(
    contig_est_ptrs <- est_ptr_from_contigs(X),
    error = function(e) {
      message("Error in est_ptr_from_contigs: ", e)
    }
  )

  tryCatch(
    sample_est_ptrs <- est_ptr_from_samples(X),
    error = function(e) {
      message("Error in est_ptr_from_samples: ", e)
    }
  )

  est_ptrs <- combine_ests(contig_est_ptrs, sample_est_ptrs)

  est_ptrs
}
