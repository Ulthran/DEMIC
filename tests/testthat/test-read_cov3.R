test_that("read_cov3 loads cov3 files as tibble", {
  tmp <- tempfile(fileext = ".cov3")
  write.csv(max_bin_001, tmp, row.names = FALSE, quote = FALSE)
  x <- read_cov3(tmp)
  expect_s3_class(x, "tbl_df")
  expect_equal(names(x), names(max_bin_001))
  expect_equal(as.data.frame(x), max_bin_001)
})
