

test_that("get_alpha_n correct calculation", {
  expect_equal(get_alpha_n(1), 0.05)
  expect_equal(get_alpha_n(1, 0.1), 0.1)
  expect_equal(get_alpha_n(100), 0.0005128014)
})



test_that("get_alpha_n wrong parameters", {
  expect_error(get_alpha_n(-1), "Sample size n not positive")
  expect_error(get_alpha_n(0), "Sample size n not positive")
  expect_error(get_alpha_n(1, 1), "alpha not beteen 0 and 1")
  expect_error(get_alpha_n(1, 2), "alpha not beteen 0 and 1")
  expect_error(get_alpha_n(1, 0), "alpha not beteen 0 and 1")
  expect_error(get_alpha_n(1, -1), "alpha not beteen 0 and 1")
})


test_that("get_robust_estimates good parameters", {
  set.seed(5)
  x <- rnorm(100)
  estimator <- get_robust_estimates(x, "norm")
  expect_equal(abs(estimator$location - 0) < 0.1, TRUE)
  expect_equal(abs(estimator$scale - 1) < 0.05, TRUE)
  x <- rlogis(100)
  estimator <- get_robust_estimates(x, "logis")
  expect_equal(abs(estimator$location - 0) < 0.4, TRUE)
  expect_equal(abs(estimator$scale - 1) < 0.2, TRUE)
  x <- rcauchy(100)
  estimator <- get_robust_estimates(x, "cauchy")
  expect_equal(abs(estimator$location - 0) < 0.2, TRUE)
  expect_equal(abs(estimator$scale - 1) < 0.1, TRUE)
  x <- rlnorm(100)
  estimator <- get_robust_estimates(x, "lnorm")
  expect_equal(abs(estimator$location - 0) < 0.15, TRUE)
  expect_equal(abs(estimator$scale - 1) < 0.15, TRUE)
})



test_that("get_robust_estimates wrong parameters", {
  x <- rnorm(100)
  expect_error(get_robust_estimates(x, "??"), "given distribution is supported yet")
})


test_that("get_bn wrong parameters", {
  expect_error(get_bn(100, alternative = "??"), "Choose alternative from greater/less/two.sided")
  expect_error(get_bn(100, distribution = "??"), "given distribution is supported yet")
  expect_error(get_bn(-1), "Sample size is not sufficient")
  expect_error(get_bn(0), "Sample size is not sufficient")
  expect_error(get_bn(2), "Sample size is not sufficient")
})



test_that("get_bn good parameters", {
  expect_equal(abs(get_bn(100) - 2.326348) < 0.01, TRUE)
  expect_equal(abs(get_bn(100, alternative = "less") - 2.326348) < 0.01, TRUE)
  expect_equal(abs(get_bn(100, alternative = "two.sided") - 2.575829) < 0.01, TRUE)
  expect_equal(abs(get_bn(100, distribution = "logis") - 4.59512) < 0.01, TRUE)
  expect_equal(abs(get_bn(100, distribution = "cauchy") - 31.83099) < 0.01, TRUE)
  expect_equal(abs(get_bn(100, distribution = "lnorm") - 10.24047) < 0.01, TRUE)
  expect_equal(abs(get_bn(100, distribution = "laplace") - 3.912023) < 0.01, TRUE)
  expect_equal(abs(get_bn(100, distribution = "gumbel") - 4.600149) < 0.01, TRUE)
})


test_that("get_an wrong parameters", {
  expect_error(get_an(100, alternative = "??"), "Choose alternative from greater/less/two.sided")
  expect_error(get_an(100, distribution = "??"), "given distribution is supported yet")
  expect_error(get_an(-1), "Sample size is not sufficient")
  expect_error(get_an(0), "Sample size is not sufficient")
  expect_error(get_an(2), "Sample size is not sufficient")
})


test_that("get_an good parameters", {
  expect_equal(abs(get_an(100) - 0.4298583) < 0.01, TRUE)
  expect_equal(abs(get_an(100, alternative = "less") - 0.4298583) < 0.01, TRUE)
  expect_equal(abs(get_an(100, alternative = "two.sided") - 0.3882245) < 0.01, TRUE)
  expect_equal(abs(get_an(100, distribution = "logis") - 1.010101) < 0.01, TRUE)
  expect_equal(abs(get_an(100, distribution = "cauchy") - 31.83099) < 0.01, TRUE)
  expect_equal(abs(get_an(100, distribution = "lnorm") - 3.84227) < 0.01, TRUE)
  expect_equal(abs(get_an(100, distribution = "laplace") - 1) < 0.01, TRUE)
  expect_equal(abs(get_an(100, distribution = "gumbel") - 1.005042) < 0.01, TRUE)
})


test_that("bp_statistic wrong parameters", {
  x <- rnorm(100)
  expect_error(bp_statistic(x, alternative = "??"), "Choose alternative from greater/less/two.sided")
  expect_error(bp_statistic(x, distribution = "??"), "given distribution is supported yet")
  expect_error(bp_statistic(x, s = 0), "Statistics limit s is not valid")
  expect_error(bp_statistic(x, s = -1), "Statistics limit s is not valid")
  expect_error(bp_statistic(x, s = 60), "Statistics limit s is not valid")

})


test_that("bp_test wrong parameters", {
  x <- rnorm(100)
  expect_error(bp_test(x, alternative = "??"), "Choose alternative from greater/less/two.sided")
  expect_error(bp_test(x, distribution = "??"), "given distribution is supported yet")

})



