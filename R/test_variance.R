#' Test the variance in PTR estimates for a given dataset
#'
#' @param X cov3 dataframe
#' @param iterations number of iterations to run
#'
#' @return named list of stats on PTR estimates
#' \itemize{
#'  \item all_sd: standard deviation of PTR estimates from all method
#'  \item all_mean: mean of PTR estimates from all method
#'  \item contigs_sd: standard deviation of PTR estimates from contigs method
#'  \item contigs_mean: mean of PTR estimates from contigs method
#'  \item samples_sd: standard deviation of PTR estimates from samples method
#'  \item samples_mean: mean of PTR estimates from samples method
#' }
#'
#' @examples
#' stats <- test_variance(ContigCluster1)
#' stats
#'
#' @export
test_variance <- function(X, iterations = 30) {
  ptrs <- list()

  for (i in 1:iterations) {
    print(paste("Iteration", i))
    ptrs[[i]] <- est_ptr(X)
  }
  print(ptrs)

  all_est_ptrs <- sapply(ptrs, function(x) x$all_ptr$est_ptr)
  contigs_est_ptrs <- sapply(ptrs, function(x) x$contigs_ptr$est_ptr)
  samples_est_ptrs <- sapply(ptrs, function(x) x$samples_ptr$est_ptr)

  list(
    all_sd = sd(all_est_ptrs),
    all_mean = mean(all_est_ptrs),
    contigs_sd = apply(contigs_est_ptrs, 1, sd),
    contigs_mean = apply(contigs_est_ptrs, 1, mean),
    samples_sd = apply(samples_est_ptrs, 1, sd),
    samples_mean = apply(samples_est_ptrs, 1, mean)
  )
}
