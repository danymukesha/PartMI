test_that("pmi returns non-negative values", {
  set.seed(42)
  x <- rnorm(100)
  y <- rnorm(100)
  z <- rnorm(100)
  expect_gte(pmi(x, y, z, method = "binning", bins = 5), 0)
  expect_gte(pmi(x, y, z, method = "knn", k = 3), 0)
})

test_that("pmi is near zero when X and Y are independent given Z", {
  set.seed(42)
  n <- 200
  z <- rnorm(n)
  x <- rnorm(n)
  y <- rnorm(n)
  pm <- pmi(x, y, z, method = "binning", bins = 5)
  expect_lt(pm, 0.2)
})

test_that("pmi reduces indirect association vs mi", {
  set.seed(42)
  n <- 200
  z <- rnorm(n)
  x <- 0.7 * z + rnorm(n, 0, sqrt(1 - 0.49))
  y <- 0.7 * z + rnorm(n, 0, sqrt(1 - 0.49))
  # x and y both depend on z, but NOT directly on each other
  mi_xy <- mi(x, y, method = "binning", bins = 5)
  pmi_xy <- pmi(x, y, z, method = "binning", bins = 5)
  expect_gt(mi_xy, pmi_xy)
})

test_that("pmi with matrix Z", {
  set.seed(42)
  n <- 100
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  x <- z1 + 0.5 * z2 + rnorm(n, 0, 0.3)
  y <- 0.5 * z1 + z2 + rnorm(n, 0, 0.3)
  z_mat <- cbind(z1, z2)
  pm <- pmi(x, y, z_mat, method = "binning", bins = 4)
  expect_type(pm, "double")
  expect_gte(pm, 0)
  pm_knn <- pmi(x, y, z_mat, method = "knn", k = 5)
  expect_type(pm_knn, "double")
  expect_gte(pm_knn, 0)
})

test_that("pmi handles NAs", {
  set.seed(42)
  x <- rnorm(100)
  y <- rnorm(100)
  z <- rnorm(100)
  x[1] <- NA
  expect_silent(pmi(x, y, z))
  expect_type(pmi(x, y, z), "double")
})

test_that("pmi with single bin sizes", {
  set.seed(42)
  x <- rnorm(50)
  y <- x + rnorm(50, 0, 0.3)
  z <- rnorm(50)
  expect_silent(pmi(x, y, z, bins = 3))
  expect_type(pmi(x, y, z, bins = 3), "double")
})
