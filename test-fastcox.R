# The fast binary-covariate risk-set core replaces coxph's stacked expansion.
# Everything downstream depends on it, so it is checked against coxph directly
# rather than trusted.

test_that("risk core reproduces coxph(ties='breslow') for a fixed covariate", {
  set.seed(11)
  for (rep in 1:5) {
    d <- simulate_trial(1200, 1200, h0_matisse, ve_shapes$gradual, 180,
                        theta = 0, ltfu_prob = 0.10, accrual_days = 60,
                        scale = 8)
    skip_if(sum(d$status) < 20)
    ref <- coxph(Surv(time, status) ~ arm, data = d, ties = "breslow")

    core <- make_risk_core(d)
    b <- 0
    for (i in 1:60) {
      p <- .core_p(core, b)
      U <- sum(core$z_ev - p); I <- sum(p * (1 - p))
      step <- U / I; b <- b + step
      if (abs(step) < 1e-12) break
    }
    p <- .core_p(core, b); se <- 1 / sqrt(sum(p * (1 - p)))

    expect_equal(unname(b), unname(coef(ref)), tolerance = 1e-6)
    expect_equal(unname(se), unname(sqrt(vcov(ref))[1, 1]), tolerance = 1e-6)
  }
})

test_that("Cox-TDC matches coxph with an explicit tt() log-time interaction", {
  set.seed(12)
  d <- simulate_trial(1500, 1500, h0_matisse, ve_shapes$rapid, 180,
                      theta = 0, ltfu_prob = 0.05, accrual_days = 90, scale = 10)
  skip_if(sum(d$status) < 40)
  ref <- coxph(Surv(time, status) ~ arm + tt(arm), data = d, ties = "breslow",
               tt = function(x, t, ...) x * log(pmax(t, 1e-8)))
  got <- est_cox_tdc(d, seq(15, 180, by = 15))

  expect_true(attr(got, "converged"))
  expect_equal(unname(attr(got, "diagnostics")$beta), unname(coef(ref)),
               tolerance = 1e-6)
  expect_equal(unname(attr(got, "diagnostics")$se_beta),
               unname(sqrt(diag(vcov(ref)))), tolerance = 1e-6)
})

test_that("the core counts the risk set correctly against a brute-force count", {
  set.seed(13)
  d <- simulate_trial(200, 200, h0_matisse, ve_shapes$constant, 180,
                      theta = 0, ltfu_prob = 0.15, accrual_days = 30, scale = 30)
  core <- make_risk_core(d)
  skip_if(core$n_ev < 5)
  for (m in seq_len(min(25, core$n_ev))) {
    tj <- core$t_ev[m]
    expect_equal(core$n0[m], sum(d$time >= tj & d$arm == 0L))
    expect_equal(core$n1[m], sum(d$time >= tj & d$arm == 1L))
  }
})
