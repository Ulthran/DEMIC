# Test case 1: Test with two vectors of equal length
test_that("lm_column returns coefficients and p-value for linear regression", {
  # Create two vectors of equal length
  x <- c(1, 2, 3, 4, 5)
  y <- c(2, 4, 6, 8, 10)

  # Call the lm_column function
  result <- lm_column(x, y)
  expect_type(result, "double")
  expect_equal(result[["y"]], 0.5)
  expect_lt(as.numeric(result[2]), 1e-6)
})
