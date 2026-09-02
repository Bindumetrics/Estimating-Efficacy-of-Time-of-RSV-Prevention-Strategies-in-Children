test_that("piecewise-constant integral and its inverse are mutually consistent", {
  pc <- h0_matisse
  t <- c(0, 1, 45, 90, 90.0001, 175, 300, 360)
  H <- pc_integral(pc, t)
  expect_true(all(diff(H) >= 0))
  expect_equal(H[1], 0)
  # inverse recovers the original times
  expect_equal(pc_integral_inverse(pc, H), t, tolerance = 1e-9)
  # Integral matches a fine-grid numeric quadrature. The tolerance is set by the
  # REFERENCE, not by pc_integral: a midpoint rule on an evenly spaced grid
  # cannot place the hazard knots exactly on cell edges, so cells straddling a
  # knot are mis-integrated by ~cell_width * delta_h, an error floor of ~1e-7
  # here. pc_integral itself is exact for a piecewise-constant integrand.
  num <- sapply(t, function(u) {
    if (u == 0) return(0)
    g <- seq(0, u, length.out = 200001)
    mean(pc_eval(pc, head(g, -1) + diff(g)[1] / 2)) * u
  })
  expect_equal(H, num, tolerance = 1e-5)
})

test_that("cumhaz_treated matches the closed form when VE is piecewise-constant", {
  brks <- c(100, 135, 170); lvls <- c(0.80, 0.65, 0.50, 0.35)
  fn <- ve_piecewise(brks, lvls)
  t <- seq(5, 360, by = 5)
  num <- cumhaz_treated(t, h0_matisse, fn)
  exa <- cumhaz_treated_exact(t, h0_matisse, brks, lvls)
  expect_equal(num, exa, tolerance = 1e-10)
})

test_that("VE_h equals VE_id when theta = 0 and differs when theta > 0", {
  t <- seq(15, 360, by = 15)
  fn <- ve_shapes$constant
  expect_equal(ve_h_true(t, h0_matisse, fn, theta = 0), fn(t))
  vh <- ve_h_true(t, h0_matisse, fn, theta = 1.0)
  expect_true(all(vh < fn(t)))            # depletion pushes VE_h below VE_id
  expect_true(all(diff(vh) < 0))          # and it wanes monotonically
})

test_that("calibration hits its target event count in expectation", {
  set.seed(7)
  for (tg in c(70, 200)) {
    sc <- calibrate_scale(tg, 3585, 3563, h0_matisse, ve_shapes$gradual, 180,
                          0, 0.04, 120)
    ev <- replicate(150, sum(simulate_trial(3585, 3563, h0_matisse,
                                            ve_shapes$gradual, 180, 0, 0.04, 120,
                                            scale = sc)$status))
    expect_lt(abs(mean(ev) - tg) / (sd(ev) / sqrt(150)), 3.5)
  }
})

test_that("realised loss to follow-up matches the nominal rate", {
  set.seed(8)
  for (p in c(0.04, 0.20)) {
    d <- simulate_trial(20000, 20000, h0_matisse, ve_shapes$gradual, 180,
                        0, ltfu_prob = p, accrual_days = 0, scale = 1)
    got <- mean(is.finite(d$t_cens) & d$t_cens < 180)
    expect_equal(got, p, tolerance = 0.01)
  }
})

test_that("administrative censoring actually censors (regression)", {
  # The original implementation put the data cutoff at accrual + horizon, which
  # makes min(horizon, A + horizon - entry) identically horizon -- the mechanism
  # censored nobody while appearing to be switched on.
  set.seed(21)
  for (f in c(0.85, 0.60)) {
    d <- simulate_trial(20000, 20000, h0_matisse, ve_shapes$gradual, 180,
                        0, ltfu_prob = 0, accrual_days = 120, scale = 1,
                        admin_complete_frac = f)
    complete <- mean(d$t_cens >= 180 - 1e-9)
    expect_equal(complete, f, tolerance = 0.02)
    # and the censored ones must be spread, not stacked at a single point
    early <- d$t_cens[d$t_cens < 180 - 1e-9]
    expect_gt(length(unique(round(early, 1))), 100)
    expect_gt(min(early), 180 - (1 - f) * 120 - 1)
  }
})

test_that("accrual on and accrual off give materially different exit patterns", {
  set.seed(22)
  a <- simulate_trial(20000, 20000, h0_matisse, ve_shapes$gradual, 180, 0, 0, 0, scale = 1)
  b <- simulate_trial(20000, 20000, h0_matisse, ve_shapes$gradual, 180, 0, 0, 120, scale = 1)
  pile_a <- mean(a$time >= 180 - 1e-9)
  pile_b <- mean(b$time >= 180 - 1e-9)
  expect_gt(pile_a - pile_b, 0.10)
})

test_that("administrative censoring alone piles every non-case at the horizon", {
  set.seed(9)
  d <- simulate_trial(4000, 4000, h0_matisse, ve_shapes$gradual, 180,
                      0, ltfu_prob = 0, accrual_days = 0, scale = 1)
  non <- d[d$status == 0L, ]
  expect_true(all(abs(non$time - 180) < 1e-9))
})

test_that("ve_logt is affine in log t on the eta scale but ve_logquad is not", {
  t <- c(10, 30, 90, 180)
  eta_t <- log(1 - ve_logt(log(0.3), 0.1)(t))
  expect_equal(eta_t, log(0.3) + 0.1 * log(t))
  # ve_logquad with c = 0 uses log(1 + t/30), which is NOT affine in log t
  eta_q <- log(1 - ve_logquad(log(0.3), 0.1, 0)(t))
  fitq <- lm(eta_q ~ log(t))
  expect_gt(max(abs(residuals(fitq))), 1e-3)
})
