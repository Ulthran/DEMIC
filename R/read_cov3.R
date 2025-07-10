#' Read a PyCov3 coverage report
#'
#' Convenient helper around readr::read_csv that converts sample and contig
#' columns to factors.
#'
#' @param file Path to a .cov3 file
#' @param ... Additional arguments passed to readr::read_csv
#'
#' @return A tibble with columns `log_cov`, `GC_content`, `sample`, `contig`,
#'   and `length`.
#'
#' @examples
#' tmp <- tempfile(fileext = ".cov3")
#' write.csv(max_bin_001, tmp, row.names = FALSE, quote = FALSE)
#' read_cov3(tmp)
#'
#' @export
read_cov3 <- function(file, ...) {
  readr::read_csv(
    file,
    col_types = readr::cols(
      log_cov = readr::col_double(),
      GC_content = readr::col_double(),
      sample = readr::col_character(),
      contig = readr::col_character(),
      length = readr::col_integer()
    ),
    ...
  ) |>
    dplyr::mutate(
      sample = factor(sample),
      contig = factor(contig)
    )
}
