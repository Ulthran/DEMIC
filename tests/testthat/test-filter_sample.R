test_that("filter_sample selects samples above thresholds", {
  Z <- data.frame(
    sample = c(rep("A", 3), rep("B", 2)),
    contig = c("c1", "c2", "c3", "c1", "c2"),
    correctY = c(5, 5, 5, 1, 1)
  )
  res <- filter_sample(Z, avg_cutoff = 2, cutoff_ratio = 0.6)
  expect_equal(res, "A")
})

test_that("filter_sample returns all samples when criteria met", {
  Z <- data.frame(
    sample = rep(c("A", "B"), each = 3),
    contig = rep(c("c1", "c2", "c3"), times = 2),
    correctY = c(3, 3, 3, 2, 2, 2)
  )
  res <- filter_sample(Z, avg_cutoff = 1, cutoff_ratio = 0.5)
  expect_setequal(res, c("A", "B"))
})
