EO1 <- structure(list(sample = c("Sample1", "Sample2", "Sample3"), est_PTR = c(1.58820988778219, 2.29699398418901, 1.84056370898365), coefficient = c(0.409285576940177, 0.73580397477001, 0.539312123134913), pValue = c(2.61567654416603e-29, 2.05157018929424e-66, 2.43569174933219e-35), cor = c(0.767467786321584, 0.93543307519703, 0.819665386998375), correctY = c(1.4891508185809, 1.59052462645664, 1.42535597156592)), row.names = c(NA, 3L), class = "data.frame")
EO2 <- structure(list(sample = c("Sample1", "Sample2", "Sample3"), est_PTR = c(2.21303380999812, 1.75222762168833, 2.25436646747005), coefficient = c(0.626249889762292, 0.442604625244878, 0.64168504060637), pValue = c(2.57415826813133e-48, 1.28302075677619e-31, 6.40172677352955e-52), cor = c(0.976659750095046, 0.928251509604572, 0.980684877999788), correctY = c(0.309937783089312, 0.377917440682641, 0.304441148907541)), row.names = c(NA, 3L), class = "data.frame")

test_that("DEMIC main function produces correct output on test set 2 cluster 1", {
  X <- read.csv("data/ContigCluster1.cov3", header = FALSE, stringsAsFactors = TRUE)
  colnames(X) <- c("log_cov", "GC_content", "sample", "contig", "length")
  O <- est_PTR(X)

  expect_equal(O, EO1)
})

test_that("DEMIC main function produces correct output on test set 2 cluster 2", {
  X <- read.csv("data/ContigCluster2.cov3", header = FALSE, stringsAsFactors = TRUE)
  colnames(X) <- c("log_cov", "GC_content", "sample", "contig", "length")
  O <- est_PTR(X)

  expect_equal(O, EO2)
})
