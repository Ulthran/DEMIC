test_that("DEMIC produces correct PTRs on generated inputs 001", {
  X <- read.csv("data/max_bin.001.cov3", header = FALSE, stringsAsFactors = TRUE)
  colnames(X) <- c("logCov", "GC", "sample", "contig", "length")
  O <- estPTR(X)

  expect_equal(O$estPTR, c(2, 3, 4), tolerance = 2)
  print(O)
})

test_that("DEMIC produces correct PTRs on generated inputs 001", {
  X <- read.csv("data/max_bin.002.cov3", header = FALSE, stringsAsFactors = TRUE)
  colnames(X) <- c("logCov", "GC", "sample", "contig", "length")
  O <- estPTR(X)

  expect_equal(O$estPTR, c(2, 3, 4), tolerance = 0.1)
  print(O)
})

test_that("DEMIC produces correct PTRs on generated inputs 001", {
  X <- read.csv("data/max_bin.003.cov3", header = FALSE, stringsAsFactors = TRUE)
  colnames(X) <- c("logCov", "GC", "sample", "contig", "length")
  O <- estPTR(X)

  expect_equal(O$estPTR, c(2, 3, 4), tolerance = 0.1)
  print(O)
})
