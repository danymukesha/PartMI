test_that("pmi_network returns expected structure", {
  set.seed(42)
  n <- 50
  g <- rnorm(n)
  x <- g + rnorm(n, 0, 0.5)
  y <- g + rnorm(n, 0, 0.5)
  z <- rnorm(n)
  dat <- data.frame(X = x, Y = y, Z = z, W = rnorm(n))
  net <- pmi_network(dat, n_permutations = 10, threshold = 0.2, verbose = FALSE)
  expect_s3_class(net, "pmi_network")
  expect_true(is.matrix(net$adjacency))
  expect_true(is.matrix(net$pmi_matrix))
  expect_true(is.matrix(net$pvalue_matrix))
  expect_true(is.matrix(net$pvalue_adjusted))
})

test_that("pmi_network adjacency matrix is symmetric", {
  set.seed(42)
  dat <- data.frame(
    A = rnorm(30),
    B = rnorm(30),
    C = rnorm(30),
    D = rnorm(30)
  )
  net <- pmi_network(dat, n_permutations = 5, threshold = 0.3, verbose = FALSE)
  expect_equal(net$adjacency, t(net$adjacency))
})

test_that("pmi_network handles edge case with 2 variables", {
  set.seed(42)
  dat <- data.frame(X = rnorm(20), Y = rnorm(20))
  net <- pmi_network(dat, n_permutations = 5, verbose = FALSE)
  expect_s3_class(net, "pmi_network")
  expect_equal(dim(net$adjacency), c(2, 2))
})

test_that("pmi_network removes zero-variance columns", {
  dat <- data.frame(A = rnorm(20), B = rep(1, 20), C = rnorm(20))
  expect_warning(
    net <- pmi_network(dat, n_permutations = 5, verbose = FALSE),
    NA
  )
  expect_equal(ncol(net$adjacency), 2)
})

test_that("print method works", {
  set.seed(42)
  dat <- data.frame(A = rnorm(20), B = rnorm(20), C = rnorm(20))
  net <- pmi_network(dat, n_permutations = 5, threshold = 0.3, verbose = FALSE)
  expect_output(print(net), "PMI Network")
})
