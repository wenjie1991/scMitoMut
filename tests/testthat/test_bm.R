## simulate data
set.seed(123)
n <- rnorm(1000, 1000, 5) %>% round(0)
theta1 <- 0.9
theta2 <- 0.1
p1 <- 0.99
p2 <- 0.8
p <- sample(c(p1, p2), 1000, replace = TRUE, prob = c(theta1, theta2)) %>% as.numeric()
x <- sapply(seq_along(n), function(i) {
  rbinom(1, n[i], p[i])
})
init_theta <- matrix(rep(0.5, 4), 2, 2)
init_lambda <- c(0.2, 0.8)

m <- data.matrix(data.frame(x, n - x))

smm_out <- em_bm(x, n, p1 = 0.90, p2 = 0.90 / 1.1, theta1 = 0.2, max_iter = 1000, tol = 1e-6)
smm_out$p2

test_that("em_bm", {
  expect_lt(abs(smm_out$theta1 - theta1) / theta1, 0.01)
  expect_lt(abs(smm_out$p1 - p1) / p1, 0.01)
  expect_lt(abs(smm_out$p2 - p2) / p2, 0.01)
})
