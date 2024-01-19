#' Estimate PTRs using all input data as well as using subsets of contigs and samples
#'
#' @param X dataframe with coverage matrix
#' (column names: "log_cov", "GC_content", "sample", "contig", "length")
#' @return named list with results from all three methods
#' all_ptr dataframe with the estimated PTRs on success, null otherwise
#' \itemize{
#'   \item est_ptr: estimated PTR values
#'   \item coefficient: coefficient of linear regression
#'   \item pValue: p-value of linear regression
#'   \item cor: correlation coefficient
#'   \item correctY: corrected coverage
#' }
#' contigs_ptr dataframe with the estimated PTRs on success, null otherwise
#' \itemize{
#'   \item est_ptr: estimated PTR values
#'   \item coefficient: coefficient of linear regression
#'   \item pValue: p-value of linear regression
#'   \item cor: correlation coefficient
#'   \item correctY: corrected coverage
#' }
#' samples_ptr dataframe with the estimated PTRs on success, null otherwise
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
  tryCatch(
    all_est_ptrs <- est_ptr_from_all(X),
    error = function(e) {
      all_est_ptrs <- NULL
      message("Error in est_ptr_from_all: ", e)
    }
  )

  tryCatch(
    contig_est_ptrs <- est_ptr_on(X, "contig"),
    error = function(e) {
      contig_est_ptrs <- NULL
      message("Error in est_ptr_on_contigs: ", e)
    }
  )

  tryCatch(
    sample_est_ptrs <- est_ptr_on(X, "sample"),
    error = function(e) {
      sample_est_ptrs <- NULL
      message("Error in est_ptr_on_samples: ", e)
    }
  )

  list(all_ptr = all_est_ptrs, contigs_ptr = contig_est_ptrs, samples_ptr = sample_est_ptrs)
}
