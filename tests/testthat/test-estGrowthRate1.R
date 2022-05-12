test_that("DEMIC main function produces out.eptr", {
  #tempdir <- normalizePath(tempdir())
  #estGrowthRate("data/all_final_contigs.cov3", paste(tempdir, "output", sep="\\"))
  estGrowthRate("data/all_final_contigs.cov3", "data/output")

  O <- read.csv(file.path(output, "out.eptr"))
  print(O)
  expect_equal("Working", "Working")
})
