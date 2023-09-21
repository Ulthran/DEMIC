test_that("DEMIC produces correct PTRs on generated inputs 001", {
  X <- read.csv("data/max_bin.001.cov3", header = FALSE, stringsAsFactors = TRUE)
  colnames(X) <- c("log_cov", "GC_content", "sample", "contig", "length")
  O <- est_ptr(X)

  expect_equal(O$est_ptr, c(2, 3, 4), tolerance = 2)
  print(O)
})

test_that("DEMIC produces correct PTRs on generated inputs 002", {
  X <- read.csv("data/max_bin.002.cov3", header = FALSE, stringsAsFactors = TRUE)
  colnames(X) <- c("log_cov", "GC_content", "sample", "contig", "length")
  O <- est_ptr(X)

  expect_equal(O$est_ptr, c(2, 3, 4), tolerance = 0.1)
  print(O)
})

test_that("DEMIC produces correct PTRs on generated inputs 003", {
  X <- read.csv("data/max_bin.003.cov3", header = FALSE, stringsAsFactors = TRUE)
  colnames(X) <- c("log_cov", "GC_content", "sample", "contig", "length")
  O <- est_ptr(X)

  expect_equal(O$est_ptr, c(2, 3, 4), tolerance = 0.1)
  print(O)
})

EO1 <- structure(list(sample = c("Sample1", "Sample2", "Sample3"), est_ptr = c(1.58820988778219, 2.29699398418901, 1.84056370898365), coefficient = c(0.409285576940177, 0.73580397477001, 0.539312123134913), pValue = c(2.61567654416603e-29, 2.05157018929424e-66, 2.43569174933219e-35), cor = c(0.767467786321584, 0.93543307519703, 0.819665386998375), correctY = c(1.4891508185809, 1.59052462645664, 1.42535597156592)), row.names = c(NA, 3L), class = "data.frame")
EO2 <- structure(list(sample = c("Sample1", "Sample2", "Sample3"), est_ptr = c(2.21303380999812, 1.75222762168833, 2.25436646747005), coefficient = c(0.626249889762292, 0.442604625244878, 0.64168504060637), pValue = c(2.57415826813133e-48, 1.28302075677619e-31, 6.40172677352955e-52), cor = c(0.976659750095046, 0.928251509604572, 0.980684877999788), correctY = c(0.309937783089312, 0.377917440682641, 0.304441148907541)), row.names = c(NA, 3L), class = "data.frame")

test_that("DEMIC main function produces correct output on sourceforge data cluster 1", {
  X <- read.csv("data/ContigCluster1.cov3", header = FALSE, stringsAsFactors = TRUE)
  colnames(X) <- c("log_cov", "GC_content", "sample", "contig", "length")
  O <- est_ptr(X)

  expect_equal(O, EO1)
})

test_that("DEMIC main function produces correct output on sourceforge data cluster 2", {
  X <- read.csv("data/ContigCluster2.cov3", header = FALSE, stringsAsFactors = TRUE)
  colnames(X) <- c("log_cov", "GC_content", "sample", "contig", "length")
  O <- est_ptr(X)

  expect_equal(O, EO2)
})
