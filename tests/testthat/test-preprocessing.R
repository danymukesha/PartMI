test_that("normalize_data zscore produces mean ~0, sd ~1", {
  set.seed(42)
  dat <- matrix(rnorm(100, mean = 5, sd = 2), nrow = 20, ncol = 5)
  nz <- normalize_data(dat, "zscore")
  expect_equal(colMeans(nz), rep(0, 5), tolerance = 1e-10)
  expect_equal(apply(nz, 2, sd), rep(1, 5), tolerance = 1e-10)
})

test_that("normalize_data minmax scales to [0,1]", {
  set.seed(42)
  dat <- matrix(runif(100, -10, 50), nrow = 20, ncol = 5)
  nm <- normalize_data(dat, "minmax")
  expect_true(all(nm >= 0 & nm <= 1))
})

test_that("normalize_data rank returns approximately normal values", {
  set.seed(42)
  dat <- matrix(rexp(100, rate = 1), nrow = 20, ncol = 5)
  nr <- normalize_data(dat, "rank")
  expect_true(abs(mean(nr)) < 0.5)
  expect_true(abs(mean(apply(nr, 2, sd)) - 1) < 0.3)
})

test_that("discretize_data returns integer matrix", {
  set.seed(42)
  dat <- matrix(rnorm(100), nrow = 20, ncol = 5)
  d <- discretize_data(dat, bins = 5)
  expect_type(d, "integer")
  expect_true(all(d >= 1 & d <= 5))
})

test_that("discretize_data quantile and width methods", {
  set.seed(42)
  dat <- matrix(rnorm(100), nrow = 20, ncol = 5)
  dq <- discretize_data(dat, method = "quantile")
  dw <- discretize_data(dat, method = "width")
  expect_type(dq, "integer")
  expect_type(dw, "integer")
  expect_equal(dim(dq), dim(dw))
})

test_that("normalize_data handles zero-variance columns", {
  dat <- cbind(A = rnorm(10), B = rep(5, 10))
  nz <- normalize_data(dat, "zscore")
  expect_true(all(is.finite(nz[, 1])))
})
