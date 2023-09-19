test_that("DEMIC produces correct PTRs on generated inputs 001", {
  X <- read.csv("data/max_bin.001.cov3", header = FALSE, stringsAsFactors = TRUE)
  colnames(X) <- c("log_cov", "GC_content", "sample", "contig", "length")
  O <- est_PTR(X)

  expect_equal(O$est_PTR, c(2, 3, 4), tolerance = 2)
  print(O)
})

test_that("DEMIC produces correct PTRs on generated inputs 001", {
  X <- read.csv("data/max_bin.002.cov3", header = FALSE, stringsAsFactors = TRUE)
  colnames(X) <- c("log_cov", "GC_content", "sample", "contig", "length")
  O <- est_PTR(X)

  expect_equal(O$est_PTR, c(2, 3, 4), tolerance = 0.1)
  print(O)
})

test_that("DEMIC produces correct PTRs on generated inputs 001", {
  X <- read.csv("data/max_bin.003.cov3", header = FALSE, stringsAsFactors = TRUE)
  colnames(X) <- c("log_cov", "GC_content", "sample", "contig", "length")
  O <- est_PTR(X)

  expect_equal(O$est_PTR, c(2, 3, 4), tolerance = 0.1)
  print(O)
})
