# ---------------------------------------------------------------------------
# GATE 2 -- estimator correctness. Every check reports a NUMBER.
# Run: Rscript scripts/gate2_estimators.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(survival) })
source("R/dgp.R"); source("R/estimators.R"); source("R/bayes.R")

set.seed(20260722)
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
res <- list(); PASS <- TRUE
h0 <- h0_matisse

FREQ <- c("cox_tdc", "durham", "tian", "pooled_logistic")
call_est <- function(m, d, grid, horizon = 180) switch(m,
  cox_tdc         = est_cox_tdc(d, grid),
  durham          = est_durham(d, grid),
  tian            = est_tian(d, grid),
  pooled_logistic = est_pooled_logistic(d, grid, interval_width = 45, horizon = horizon),
  bayes           = est_bayes(d, grid, horizon = horizon, seed = 7))

# ---------------------------------------------------------------------------
say("\n=== 2.1 Monotone-transform / interval-inversion check ===")
say("  VE = 1 - exp(eta) is DECREASING in eta, so the upper eta bound must give")
say("  the LOWER VE bound. A reversed mapping inflates coverage plausibly.")
sc <- calibrate_scale(200, 3585, 3563, h0, ve_shapes$gradual, 180, 0, 0.04, 120)
d  <- simulate_trial(3585, 3563, h0, ve_shapes$gradual, 180, 0, 0.04, 120, scale = sc)
grid <- seq(15, 180, by = 15)
ok21 <- TRUE
for (m in c(FREQ, "bayes")) {
  r <- call_est(m, d, grid)
  v <- !is.na(r$ve_hat)
  ord <- all(r$lo[v] <= r$ve_hat[v] + 1e-12) && all(r$ve_hat[v] <= r$hi[v] + 1e-12)
  wid <- all(r$hi[v] - r$lo[v] > 0)
  say("  %-16s lo <= est <= hi : %-5s   all widths positive : %-5s   (%d/%d estimable)",
      m, ord, wid, sum(v), length(v))
  ok21 <- ok21 && ord && wid
}
# direct algebraic assertion on the shared mapping
tst <- ve_from_eta(eta = c(-1.2, -0.3), se = c(0.4, 0.25))
alg <- all(tst$lo < tst$ve_hat) && all(tst$ve_hat < tst$hi) &&
  isTRUE(all.equal(tst$lo, 1 - exp(c(-1.2, -0.3) + qnorm(0.975) * c(0.4, 0.25))))
say("  ve_from_eta() algebraic identity holds: %s", alg)
PASS <- PASS && ok21 && alg
say("  [%s] interval inversion correct for all five methods", if (ok21 && alg) "PASS" else "FAIL")

# ---------------------------------------------------------------------------
say("\n=== 2.2 Known-answer test: constant VE, proportional hazards ===")
say("  Under constant VE every method must recover the constant log HR, and the")
say("  Cox-TDC log-t interaction must be indistinguishable from zero.")
say("  Scored over REPLICATES, not on one dataset. On a single draw with ~4000")
say("  events the sampling SE of VE is ~0.011, so a -0.013 deviation is one")
say("  standard error of noise and says nothing about bias either way.")
ve0 <- 0.70; eta_true <- log(1 - ve0)
sc <- calibrate_scale(4000, 20000, 20000, h0, ve_shapes$constant, 180, 0, 0.04, 120)
nrep22 <- 60
bias_mat <- matrix(NA, nrep22, length(FREQ), dimnames = list(NULL, FREQ))
cox_b2 <- numeric(nrep22); cox_p <- numeric(nrep22)
for (i in seq_len(nrep22)) {
  dd <- simulate_trial(20000, 20000, h0, ve_shapes$constant, 180, 0, 0.04, 120, scale = sc)
  for (m in FREQ) bias_mat[i, m] <- mean(call_est(m, dd, grid)$ve_hat - ve0, na.rm = TRUE)
  dc <- attr(est_cox_tdc(dd, grid), "diagnostics")
  cox_b2[i] <- dc$beta[2]; cox_p[i] <- dc$wald_p_waning
}
say("  %d replicates at n = 40000 (~4000 events each; the check needs EVENTS,
  not subjects, so the sample was reduced and the hazard scale raised)", nrep22)
ok22 <- TRUE
for (m in FREQ) {
  b <- mean(bias_mat[, m]); s <- sd(bias_mat[, m]) / sqrt(nrep22)
  z <- b / s
  say("  %-16s mean bias = %+.5f   MCSE = %.5f   z = %+5.2f", m, b, s, z)
  res[[paste0("known.bias.", m)]] <- b
  res[[paste0("known.z.", m)]] <- z
  ok22 <- ok22 && abs(z) < 3
}
say("  Cox-TDC log-t interaction b2: mean = %+.5f (MCSE %.5f); truth = 0",
    mean(cox_b2), sd(cox_b2) / sqrt(nrep22))
say("  type I error of the no-waning Wald test at 5%% = %.1f%% (MCSE %.1f%%)",
    100 * mean(cox_p < 0.05), 100 * sqrt(mean(cox_p < 0.05) * (1 - mean(cox_p < 0.05)) / nrep22))
res[["known.typeI_cox_tdc"]] <- mean(cox_p < 0.05)
ok22 <- ok22 && abs(mean(cox_b2) / (sd(cox_b2) / sqrt(nrep22))) < 3
PASS <- PASS && ok22
say("  [%s] all methods recover constant VE; no spurious waning detected",
    if (ok22) "PASS" else "FAIL")

# ---------------------------------------------------------------------------
say("\n=== 2.3 Consistency: bias must shrink with n ===")
say("  Consistency is tested on a CONSTANT-VE truth, where all four methods are")
say("  exactly specified. On a curved truth, Durham (fixed df = 4) and pooled")
say("  logistic (fixed 45-day intervals) converge to a nonzero SMOOTHING bias as")
say("  n grows -- they are consistent only as their tuning refines, not as n")
say("  alone increases. Requiring their bias to vanish with n at fixed tuning")
say("  would be testing a property they do not have and do not claim.")
fn <- ve_shapes$constant
tab <- list(); tabse <- list()
for (nn in c(1e4, 1e5)) {
  n_arm <- nn / 2
  scn <- calibrate_scale(nn * 200 / 7148, n_arm, n_arm, h0, fn, 180, 0, 0.04, 120)
  nrep23 <- if (nn == 1e4) 60 else 40
  b <- sapply(seq_len(nrep23), function(i) {
    dd <- simulate_trial(n_arm, n_arm, h0, fn, 180, 0, 0.04, 120, scale = scn)
    sapply(FREQ, function(m) mean(call_est(m, dd, grid)$ve_hat - fn(grid), na.rm = TRUE))
  })
  tab[[as.character(nn)]] <- rowMeans(b, na.rm = TRUE)
  tabse[[as.character(nn)]] <- apply(b, 1, sd, na.rm = TRUE) / sqrt(nrep23)
  say("  n = %6.0e (%d reps): %s", nn, nrep23,
      paste(sprintf("%s %+0.4f(%.4f)", FREQ, rowMeans(b, na.rm = TRUE),
                    apply(b, 1, sd, na.rm = TRUE) / sqrt(nrep23)), collapse = "  "))
}
ok23 <- TRUE
say("  bias with MCSE in brackets; z = bias/MCSE at the larger n:")
for (m in FREQ) {
  r10 <- abs(tab[["10000"]][m]); r100 <- abs(tab[["1e+05"]][m])
  z <- tab[["1e+05"]][m] / tabse[["1e+05"]][m]
  say("  %-16s |bias| %.5f -> %.5f   z at n=1e5 = %+5.2f", m, r10, r100, z)
  res[[paste0("consistency.z.", m)]] <- z
  ok23 <- ok23 && abs(z) < 3
}
PASS <- PASS && ok23
say("  [%s] bias indistinguishable from zero at n = 1e5 for all four",
    if (ok23) "PASS" else "FAIL")

say("\n  Companion panel (NOT a pass/fail): irreducible smoothing bias on a")
say("  curved truth at fixed tuning, n = 1e5 -- evidence for the")
say("  'consistent-but-smoothed' tier of the specification map.")
fnc <- ve_shapes$gradual
scn <- calibrate_scale(1e5 * 200 / 7148, 5e4, 5e4, h0, fnc, 180, 0, 0.04, 120)
bs <- sapply(1:40, function(i) {
  dd <- simulate_trial(5e4, 5e4, h0, fnc, 180, 0, 0.04, 120, scale = scn)
  sapply(FREQ, function(m) mean(call_est(m, dd, grid)$ve_hat - fnc(grid), na.rm = TRUE))
})
for (m in FREQ) {
  b <- mean(bs[m, ]); s <- sd(bs[m, ]) / sqrt(ncol(bs))
  say("    %-16s bias = %+.5f (MCSE %.5f)  z = %+6.2f", m, b, s, b / s)
  res[[paste0("smoothing_bias.", m)]] <- b
}

# ---------------------------------------------------------------------------
say("\n=== 2.4 Interval calibration where a method is correctly specified ===")
say("  Cox-TDC is EXACT when log(1-VE(t)) is affine in log t. Build exactly that")
say("  truth, give it abundant data, and its pointwise coverage must be ~95%%.")
say("  NB this uses ve_logt(), the Cox-TDC family itself. An earlier version of")
say("  this check used ve_logquad(c = 0), whose basis is log(1 + t/30) -- NOT an")
say("  affine function of log t -- so Cox-TDC was misspecified and the check was")
say("  measuring misspecification bias (4.8 MCSE under-coverage), not calibration.")
fn_aff <- ve_logt(a = log(0.30), b = 0.10)
nrep <- 1000
cov_mat <- matrix(NA, nrep, length(grid))
sc_aff <- calibrate_scale(1500, 20000, 20000, h0, fn_aff, 180, 0, 0.04, 120)
truth_aff <- fn_aff(grid)
for (i in seq_len(nrep)) {
  dd <- simulate_trial(20000, 20000, h0, fn_aff, 180, 0, 0.04, 120, scale = sc_aff)
  r <- est_cox_tdc(dd, grid)
  cov_mat[i, ] <- as.numeric(r$lo <= truth_aff & truth_aff <= r$hi)
}
cvg <- colMeans(cov_mat, na.rm = TRUE)
mcse <- sqrt(cvg * (1 - cvg) / nrep)
say("     t   coverage   MCSE")
for (k in seq_along(grid)) say("  %4d     %5.1f%%   %4.1f%%", grid[k], 100 * cvg[k], 100 * mcse[k])
worst <- max(abs(cvg - 0.95))
say("  mean coverage = %.1f%%, max |coverage - 95%%| = %.1f pp (typical MCSE %.1f pp)",
    100 * mean(cvg), 100 * worst, 100 * mean(mcse))
res[["coverage.cox_tdc.max_dev"]] <- worst
res[["coverage.cox_tdc.max_z"]] <- max(abs(cvg - 0.95) / mcse)
say("  worst pointwise deviation = %.2f MCSE", max(abs(cvg - 0.95) / mcse))
ok24 <- max(abs(cvg - 0.95) / mcse) < 3.5
PASS <- PASS && ok24
say("  [%s] coverage within Monte Carlo error of nominal", if (ok24) "PASS" else "FAIL")

# ---------------------------------------------------------------------------
say("\n=== 2.5 Cross-estimator agreement (D'Agostino equivalence) ===")
say("  Two separate checks. (a) The equivalence D'Agostino actually proves is for")
say("  a TIME-CONSTANT effect with short intervals: cloglog grouped survival ->")
say("  Cox as the interval width goes to zero. (b) For the time-VARYING versions")
say("  the two are different estimators, so they are compared as a PAIRED mean")
say("  difference over replicates -- a single-draw difference is mostly noise.")
say("")
say("  (a) constant-effect equivalence, interval width -> 0:")
dd <- simulate_trial(20000, 20000, h0, ve_shapes$constant, 180, 0, 0.04, 120, scale = sc)
ref_c <- coef(coxph(Surv(time, status) ~ arm, data = dd, ties = "breslow"))
for (w in c(60, 30, 15, 5)) {
  brk <- seq(0, 180, by = w); K <- length(brk) - 1L
  rows <- lapply(seq_len(K), function(k) {
    lo <- brk[k]; hi <- brk[k + 1]; ar <- dd$time > lo
    if (!any(ar)) return(NULL)
    s <- dd[ar, ]
    data.frame(interval = factor(k, levels = seq_len(K)), arm = s$arm,
               event = as.integer(s$status == 1L & s$time <= hi),
               logexp = log(pmax(pmin(s$time, hi) - lo, 1e-6)))
  })
  pl <- do.call(rbind, rows[!vapply(rows, is.null, TRUE)])
  g <- suppressWarnings(glm(event ~ interval - 1 + arm + offset(logexp),
                            family = binomial("cloglog"), data = pl))
  say("     width %3d d (K=%2d): cloglog beta = %+.5f   coxph = %+.5f   diff = %+.6f",
      w, K, coef(g)[["arm"]], ref_c, coef(g)[["arm"]] - ref_c)
  if (w == 5) { dev5 <- abs(coef(g)[["arm"]] - ref_c); res[["agreement.constant_diff"]] <- dev5 }
}
say("")
say("  (b) paired mean difference of the time-varying versions, %d replicates:", nrep22)
difs <- matrix(NA, nrep22, length(grid))
for (i in seq_len(nrep22)) {
  dd2 <- simulate_trial(20000, 20000, h0, ve_shapes$constant, 180, 0, 0.04, 120, scale = sc)
  difs[i, ] <- est_cox_tdc(dd2, grid)$ve_hat -
    est_pooled_logistic(dd2, grid, interval_width = 45, horizon = 180)$ve_hat
}
md <- colMeans(difs, na.rm = TRUE)
sd_ <- apply(difs, 2, sd, na.rm = TRUE) / sqrt(nrep22)
say("     max |mean paired difference| = %.5f VE units (largest MCSE %.5f, max |z| = %.2f)",
    max(abs(md)), max(sd_), max(abs(md / sd_)))
say("     mean |per-replicate difference| = %.5f -- this is the sampling spread,",
    mean(abs(difs), na.rm = TRUE))
say("     which is what a single-draw comparison would have been measuring.")
res[["agreement.paired_max_z"]] <- max(abs(md / sd_))
ok25 <- dev5 < 0.01 && max(abs(md / sd_)) < 4
PASS <- PASS && ok25
say("  [%s] constant-effect equivalence holds and no systematic disagreement",
    if (ok25) "PASS" else "FAIL")

# ---------------------------------------------------------------------------
say("\n=== 2.6 Realised degrees of freedom in the Durham spline ===")
rd <- est_durham(d, grid)
dgd <- attr(rd, "diagnostics")
say("  realised spline rank = %d (must be exactly 4: ns(df=4, intercept=TRUE))",
    dgd$realised_df)
say("  Grambsch-Therneau variance scale used = %.6g", dgd$gt_var)
say("  grid points extrapolated beyond event support = %d", sum(dgd$extrapolated))
ok26 <- dgd$realised_df == 4
PASS <- PASS && ok26
say("  [%s] no accidental extra intercept", if (ok26) "PASS" else "FAIL")

# ---------------------------------------------------------------------------
say("\n=== 2.7 Edge cases ===")
mk <- function(n, ev_arm, ntimes) {
  data.frame(id = 1:n, arm = rep(0:1, length.out = n),
             time = c(seq_len(ntimes) * 10, rep(180, n - ntimes)),
             status = c(rep(1L, ntimes), rep(0L, n - ntimes)))
}
cases <- list(
  "zero events"          = transform(mk(200, 0, 0), status = 0L),
  "single event time"    = mk(200, 1, 1),
  "all events in one arm"= within(mk(200, 1, 12), { arm[status == 1L] <- 0L }),
  "heavy ties"           = within(mk(400, 1, 40), { time[status == 1L] <- 60 })
)
for (nm in names(cases)) {
  cd <- cases[[nm]]
  line <- sapply(FREQ, function(m) {
    r <- try(call_est(m, cd, grid), silent = TRUE)
    if (inherits(r, "try-error")) return("THREW")
    if (!attr(r, "converged")) return("fail(clean)")
    sprintf("%d/%d NA", sum(is.na(r$ve_hat)), length(grid))
  })
  say("  %-22s %s", nm, paste(sprintf("%s=%s", FREQ, line), collapse = "  "))
}
threw <- any(sapply(names(cases), function(nm)
  any(sapply(FREQ, function(m) inherits(try(call_est(m, cases[[nm]], grid), silent = TRUE),
                                        "try-error")))))
PASS <- PASS && !threw
say("  [%s] no estimator throws an uncaught error on a degenerate input",
    if (!threw) "PASS" else "FAIL")

# ---------------------------------------------------------------------------
say("\n=== GATE 2 SUMMARY ===")
say("  Overall: %s", if (PASS) "PASS" else "FAIL")
write.csv(data.frame(check = names(res), value = unlist(res)),
          "outputs/environment/gate2_results.csv", row.names = FALSE)
say("  Numbers written to outputs/environment/gate2_results.csv")
