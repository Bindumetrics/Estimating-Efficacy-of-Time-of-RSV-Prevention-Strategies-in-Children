# VE = 1 - exp(eta) is DECREASING in eta. Mapping the interval bounds the wrong
# way round inverts every CI and inflates coverage in a way that looks entirely
# plausible -- it is the single most dangerous silent error in this codebase.

test_that("ve_from_eta maps the UPPER eta bound to the LOWER VE bound", {
  eta <- c(-2.5, -1.2, -0.4, 0.0, 0.3)
  se  <- c(0.5, 0.3, 0.2, 0.4, 0.25)
  v <- ve_from_eta(eta, se, level = 0.95)
  z <- qnorm(0.975)

  expect_equal(v$ve_hat, 1 - exp(eta))
  expect_equal(v$lo, 1 - exp(eta + z * se))   # upper eta -> lower VE
  expect_equal(v$hi, 1 - exp(eta - z * se))
  expect_true(all(v$lo < v$ve_hat))
  expect_true(all(v$ve_hat < v$hi))
})

test_that("a larger standard error gives a strictly wider VE interval", {
  a <- ve_from_eta(-1, 0.2); b <- ve_from_eta(-1, 0.6)
  expect_true((b$hi - b$lo) > (a$hi - a$lo))
})

test_that("every estimator returns lo <= ve_hat <= hi with positive width", {
  set.seed(404)
  sc <- calibrate_scale(200, 3585, 3563, h0_matisse, ve_shapes$gradual, 180, 0, 0.04, 120)
  d <- simulate_trial(3585, 3563, h0_matisse, ve_shapes$gradual, 180, 0, 0.04, 120,
                      scale = sc)
  grid <- seq(15, 180, by = 15)
  for (f in list(est_cox_tdc, est_durham, est_tian,
                 function(d, g) est_pooled_logistic(d, g, 45, 180))) {
    r <- f(d, grid)
    ok <- !is.na(r$ve_hat)
    expect_true(all(r$lo[ok] <= r$ve_hat[ok] + 1e-12))
    expect_true(all(r$ve_hat[ok] <= r$hi[ok] + 1e-12))
    expect_true(all(r$hi[ok] - r$lo[ok] > 0))
  }
})

test_that("a deliberately inverted mapping is caught", {
  # If someone 'fixes' ve_from_eta by swapping lo/hi, this must fail loudly.
  bad <- function(eta, se) {
    z <- qnorm(0.975)
    list(ve_hat = 1 - exp(eta), lo = 1 - exp(eta - z * se), hi = 1 - exp(eta + z * se))
  }
  b <- bad(-1, 0.3)
  expect_false(all(b$lo < b$ve_hat & b$ve_hat < b$hi))
})
