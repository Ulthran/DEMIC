# test_that("DEMIC main function produces out.eptr on test set 1", {
# tempdir <- normalizePath(tempdir())
# estGrowthRate("data/all_final_contigs.cov3", paste(tempdir, "output", sep="\\"))
#  estGrowthRate("data/all_final_contigs.cov3", "data/output")

#  O <- read.table(file = file.path("data/output", "out.eptr"), sep = "\t", header = TRUE)
#  E <- read.table(file = file.path("data/expected_output", "out.eptr"), sep = "\t", header = TRUE)
#  expect_equal(O, E)
# })

test_that("DEMIC main function produces out.eptr on test set 2", {
  # tempdir <- normalizePath(tempdir())
  # estGrowthRate("data2/all_final_contigs.cov3", paste(tempdir, "output", sep="\\"))
  estGrowthRate("data2/ContigCluster1.cov3", output = "data2/output", log_level = DEBUG)

  O <- read.table(file = file.path("data2/output", "out.eptr"), sep = "\t", header = TRUE)
  O[, -1] <- round(O[, -1], 5)
  E <- read.table(file = file.path("data2/expected_output", "ContigCluster1.ptr"), sep = "\t", header = TRUE)
  E[, -1] <- round(E[, -1], 5)
  expect_equal(O, E)


  estGrowthRate("data2/ContigCluster2.cov3", output = "data2/output", log_level = DEBUG)

  O <- read.table(file = file.path("data2/output", "out.eptr"), sep = "\t", header = TRUE)
  O[, -1] <- round(O[, -1], 5)
  E <- read.table(file = file.path("data2/expected_output", "ContigCluster2.ptr"), sep = "\t", header = TRUE)
  E[, -1] <- round(E[, -1], 5)
  expect_equal(O, E)
})
