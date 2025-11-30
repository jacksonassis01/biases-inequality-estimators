test_that("compara o valor esperado do estimador TT com o valor real quando alpha1 = alpha2", {
  n = 15
  alpha = .5
  trueValue = \(n, alpha) digamma(alpha) - digamma(n * alpha) + (n-1) / (n*alpha) + log(n)
  expect_equal(compute_expected_value_TT_estimator(n, c(.25, .75), c(alpha, alpha)), trueValue(n, alpha))
  alpha = 2
  expect_equal(compute_expected_value_TT_estimator(n, c(.25, .75), c(alpha, alpha)), trueValue(n, alpha))
})

test_that("compara o valor esperado do estimador TL com o valor real quando alpha1 = alpha2", {
  n = 15
  alpha = .5
  trueValue = \(n, alpha) digamma(n * alpha) - digamma(alpha) - log(n)
  expect_equal(compute_expected_value_TL_estimator(n, c(.25, .75), c(alpha, alpha)), trueValue(n, alpha))
  alpha = 2
  expect_equal(compute_expected_value_TL_estimator(n, c(.25, .75), c(alpha, alpha)), trueValue(n, alpha))
})

test_that("compara o valor esperado do estimador A1 com o valor real quando alpha1 = alpha2", {
  n = 15
  alpha = .5
  trueValue = \(n, alpha) 1 - 1/alpha*(gamma(alpha + 1/n) / gamma(alpha))^n
  expect_equal(compute_expected_value_A1_estimator(n, c(.25, .75), c(alpha, alpha)), trueValue(n, alpha))
  alpha = 2
  expect_equal(compute_expected_value_A1_estimator(n, c(.25, .75), c(alpha, alpha)), trueValue(n, alpha))
})

test_that("compara o valor esperado do estimador AInf com o valor real quando alpha1 = alpha2", {
  n = 15
  alpha = .5
  trueValue <- \(n, alpha) {
    Q_gamma <- \(u, alpha) pgamma(u, shape = alpha, lower.tail = FALSE)
    integ   <- integrate(\(u) Q_gamma(u, alpha)^n, lower = 0, upper = Inf)
    return(1 - (integ$value / alpha))
  }

  expect_equal(compute_expected_value_AInf_estimator(n, c(.25, .75), c(alpha, alpha)), trueValue(n, alpha))
  alpha = 2
  expect_equal(compute_expected_value_AInf_estimator(n, c(.25, .75), c(alpha, alpha)), trueValue(n, alpha))
})

test_that("compara o valor esperado do estimador VMR com o valor real quando alpha1 = alpha2", {
  n=15; beta=1
  alpha = .5
  trueValue = \(n, alpha, beta) n * alpha / (beta * (n * alpha + 1))
  expect_equal(compute_expected_value_VMR_estimator(n, c(.25, .75), c(alpha, alpha), beta), trueValue(n, alpha, beta))
  alpha = 2
  expect_equal(compute_expected_value_VMR_estimator(n, c(.25, .75), c(alpha, alpha), beta), trueValue(n, alpha, beta))
})
