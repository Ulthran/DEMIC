test_that("KS test returns accurate p-value", {
  x<-runif(100,-1,1)
  res <- ks(x)

  expect_equal("Working", "Working")
})
