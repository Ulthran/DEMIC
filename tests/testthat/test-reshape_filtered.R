test_that("reshape_filtered drops contigs with missing values", {
  Z <- data.frame(
    sample = c("A", "A", "A", "B", "B"),
    contig = c("c1", "c2", "c3", "c1", "c2"),
    correctY = c(1, 2, 3, 4, 5)
  )
  res <- reshape_filtered(c("A", "B"), Z)
  expect_equal(rownames(res), c("c1", "c2"))
  expect_equal(colnames(res), c("A", "B"))
  expect_equal(as.matrix(res), matrix(c(1, 2, 4, 5), nrow = 2, byrow = FALSE,
                                     dimnames = list(c("c1", "c2"), c("A", "B"))))
})
