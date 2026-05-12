test_that("mi returns non-negative values", {
  set.seed(42)
  x <- rnorm(100)
  y <- rnorm(100)
  expect_gte(mi(x, y, method = "binning"), 0)
  expect_gte(mi(x, y, method = "knn"), 0)
})

test_that("mi is symmetric", {
  set.seed(42)
  x <- rnorm(100)
  y <- x + rnorm(100, 0, 0.5)
  expect_equal(mi(x, y, method = "binning"), mi(y, x, method = "binning"),
               tolerance = 1e-10)
  expect_equal(mi(x, y, method = "knn"), mi(y, x, method = "knn"),
               tolerance = 1e-10)
})

test_that("mi is zero for independent variables", {
  set.seed(42)
  x <- rnorm(200)
  y <- rnorm(200)
  m <- mi(x, y, method = "binning", bins = 5)
  expect_lt(m, 0.1)
})

test_that("mi is high for dependent variables", {
  set.seed(42)
  x <- rnorm(100)
  y <- x
  m <- mi(x, y, method = "binning", bins = 10)
  expect_gt(m, 0.5)
})

test_that("mi handles NAs", {
  set.seed(42)
  x <- rnorm(100)
  y <- rnorm(100)
  x[1] <- NA
  y[2] <- NA
  expect_silent(mi(x, y))
  expect_type(mi(x, y), "double")
})

test_that("mi binning with different bins", {
  set.seed(42)
  x <- rnorm(100)
  y <- x + rnorm(100, 0, 0.3)
  m3 <- mi(x, y, method = "binning", bins = 3)
  m10 <- mi(x, y, method = "binning", bins = 10)
  expect_type(m3, "double")
  expect_type(m10, "double")
})

test_that("mi knn with different k", {
  set.seed(42)
  x <- rnorm(100)
  y <- x + rnorm(100, 0, 0.3)
  m <- mi(x, y, method = "knn", k = 5)
  expect_type(m, "double")
  expect_gte(m, 0)
})
