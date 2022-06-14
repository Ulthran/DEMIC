test_that("DEMIC main function produces out.eptr", {
  #tempdir <- normalizePath(tempdir())
  #estGrowthRate("data/all_final_contigs.cov3", paste(tempdir, "output", sep="\\"))
  estGrowthRate("data/all_final_contigs.cov3", "data/output")

  O <- read.csv(file.path("data/output", "out.eptr"))
  E <- read.csv(file.path("data/expected_output", "out.eptr"))
  expect_equal(O, E)
})
