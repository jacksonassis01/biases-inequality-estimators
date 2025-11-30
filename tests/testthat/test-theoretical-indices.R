test_that("compute_Gini_G coincide com o valor exato quando alpha1 = alpha2", {
  alpha = .5
  trueValue = \(alpha) gamma(alpha + 0.5) / (sqrt(pi) * alpha * gamma(alpha))
  expect_equal(compute_Gini_G(pii = c(.75, .25), alpha = c(alpha, alpha)), trueValue(alpha))
  alpha = 2
  expect_equal(compute_Gini_G(pii = c(.75, .25), alpha = c(alpha, alpha)), trueValue(alpha))
})

test_that("compute_TheilT coincide com o valor exato quando alpha1 = alpha2", {
  alpha = .5
  trueValue = \(alpha) digamma(alpha) + 1/alpha - log(alpha)
  expect_equal(compute_TheilT(pii = c(.75, .25), alpha = c(alpha, alpha)), trueValue(alpha))
  alpha = 2
  expect_equal(compute_TheilT(pii = c(.75, .25), alpha = c(alpha, alpha)), trueValue(alpha))
})

test_that("compute_TheilL coincide com o valor exato quando alpha1 = alpha2", {
  alpha = .5
  trueValue = \(alpha) log(alpha) - digamma(alpha)
  expect_equal(compute_TheilL(pii = c(.75, .25), alpha = c(alpha, alpha)), trueValue(alpha))
  alpha = 2
  expect_equal(compute_TheilL(pii = c(.75, .25), alpha = c(alpha, alpha)), trueValue(alpha))
})

test_that("compute_Atkinson1 coincide com o valor exato quando alpha1 = alpha2", {
  alpha = .5
  trueValue = \(alpha) 1 - 1/alpha * exp(digamma(alpha))
  expect_equal(compute_Atkinson1(pii = c(.75, .25), alpha = c(alpha, alpha)), trueValue(alpha))
  alpha = 2
  expect_equal(compute_Atkinson1(pii = c(.75, .25), alpha = c(alpha, alpha)), trueValue(alpha))
})

test_that("compute_VMR retorna 1 quando alpha1 = alpha2", {
  expect_equal(compute_VMR(pii = c(.75, .25), alpha = c(.5, .5), beta = 1), 1)
  expect_equal(compute_VMR(pii = c(.75, .25), alpha = c(2, 2), beta = 1), 1)
})
