test_that("log beta binomial pmf", {
  a <- 10
  b <- 10
  x <- 10
  n <- 20
  expect_equal(log_beta_binomial_pmf(a, b, x, n), log(VGAM::dzoibetabinom.ab(x, n, a, b)))
  a <- 20
  b <- 10
  x <- 10
  n <- 30
  expect_equal(log_beta_binomial_pmf(a, b, x, n), log(VGAM::dzoibetabinom.ab(x, n, a, b)))
})


test_that("pbetabinom", {
  a <- 10
  b <- 10
  x <- 10
  n <- 20
  expect_equal(pbetabinom(x, n, a, b), VGAM::pbetabinom.ab(x, n, a, b))
  a <- 20
  b <- 10
  x <- 10
  n <- 30
  expect_equal(pbetabinom(x, n, a, b), VGAM::pbetabinom.ab(x, n, a, b))
  a <- 20
  b <- 10000
  x <- 10
  n <- 5000
  expect_equal(pbetabinom(x, n, a, b), VGAM::pbetabinom.ab(x, n, a, b))
  a <- 20
  b <- 10000
  x <- 50
  n <- 1000
  expect_equal(pbetabinom(x, n, a, b), VGAM::pbetabinom.ab(x, n, a, b))
})

test_that("log_n_choose_k", {
  expect_equal(log_n_choose_k(10, 2), log(choose(10, 2)))
  expect_equal(log_n_choose_k(10, 10), log(choose(10, 10)))
  expect_equal(log_n_choose_k(10, 7), log(choose(10, 7)))
  expect_equal(log_n_choose_k(1000, 7), log(choose(1000, 7)))
})


#######################################################################
#                         Test model fitting                          #
#######################################################################

## Simulate
set.seed(123)
n <- round(rnorm(10000, 70, sd = 5))
p <- 0.9
t <- 10
(s1 <- p * t)
(s2 <- t - p * t)
# s1 = 10
# s2 = 20
N <- length(n)
# (pi = s1/(s1 + s2))
# (sigma = s1 + s2)
x <- sapply(n, function(x) {
  VGAM::rbetabinom.ab(n = 1, size = x, shape1 = s1, shape2 = s2)
})

test_that("fit bb", {
  y <- mle_bb(x, n)
  expect_lt(abs(y$a / (y$a + y$b) - p) / p, 0.05)
  expect_lt(abs(y$a + y$b - t) / t, 0.05)
})
