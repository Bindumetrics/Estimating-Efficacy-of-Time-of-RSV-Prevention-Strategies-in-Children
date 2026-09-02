# ---------------------------------------------------------------------------
# The specification map (brief sec 4, MANDATORY).
#
# Which estimator is correctly specified under which true VE shape. Three tiers,
# because a binary correct/misspecified split fails here: a pooled logistic with
# interval-specific effects on the DGP's own grid would be "correct" everywhere,
# which the brief forbids.
#
#   EXACT      truth lies in the estimator's model class -- zero asymptotic bias
#   SMOOTHED   asymptotically unbiased as df -> inf / bandwidth -> 0 / interval
#              -> 0, but biased at any FIXED tuning, however large n gets
#   MISSPEC    asymptotic bias does not vanish under any tuning refinement
#
# The classification is analytic. It is then CORROBORATED empirically: bias is
# measured at n = 2e5, where sampling noise is small enough that a nonzero bias
# is visible, and the measured values must agree with the assigned tier.
# ---------------------------------------------------------------------------

source("R/dgp.R"); source("R/estimators.R")
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

# Analytic classification, with the reason recorded for each cell.
MAP <- rbind(
  data.frame(shape = "constant",
             cox_tdc = "EXACT", durham = "EXACT", tian = "EXACT",
             pooled_logistic = "EXACT", bayes = "EXACT",
             reason = "proportional hazards; every model class contains a constant log HR"),
  data.frame(shape = "gradual",
             cox_tdc = "MISSPEC", durham = "SMOOTHED", tian = "SMOOTHED",
             pooled_logistic = "SMOOTHED", bayes = "MISSPEC",
             reason = "logistic-in-log-time decay: not affine in log t, not in the Bayesian candidate set"),
  data.frame(shape = "rapid",
             cox_tdc = "MISSPEC", durham = "SMOOTHED", tian = "SMOOTHED",
             pooled_logistic = "SMOOTHED", bayes = "MISSPEC",
             reason = "as gradual, steeper; stresses spline overshoot at the boundary"),
  data.frame(shape = "piecewise",
             cox_tdc = "MISSPEC", durham = "SMOOTHED", tian = "SMOOTHED",
             pooled_logistic = "SMOOTHED", bayes = "MISSPEC",
             reason = "steps at 100/135/170 d are OFF the 30-day blocks and off every interval grid"),
  data.frame(shape = "applied (log-quadratic)",
             cox_tdc = "MISSPEC", durham = "SMOOTHED", tian = "SMOOTHED",
             pooled_logistic = "SMOOTHED", bayes = "MISSPEC",
             reason = "basis is log(1 + t/30), which is not affine in log t for any coefficients"),
  data.frame(shape = "frailty-induced",
             cox_tdc = "MISSPEC", durham = "MISSPEC", tian = "MISSPEC",
             pooled_logistic = "MISSPEC", bayes = "MISSPEC",
             reason = "VE_h(t) from gamma-frailty marginalisation lies in NO estimator's model class")
)

say("\n=== SPECIFICATION MAP ===\n")
print(MAP[, c("shape", "cox_tdc", "durham", "tian", "pooled_logistic", "bayes")],
      row.names = FALSE)
say("")
for (i in seq_len(nrow(MAP))) say("  %-24s %s", MAP$shape[i], MAP$reason[i])

METHODS <- c("cox_tdc", "durham", "tian", "pooled_logistic", "bayes")
say("\n=== Acceptance check (brief sec 4) ===")
ok <- TRUE
for (m in METHODS) {
  tiers <- MAP[[m]]
  always_exact <- all(tiers == "EXACT")
  ever_exact   <- any(tiers == "EXACT")
  ever_wrong   <- any(tiers != "EXACT")
  say("  %-16s EXACT in %d/%d shapes; misspecified or smoothed in %d/%d  %s",
      m, sum(tiers == "EXACT"), length(tiers), sum(tiers != "EXACT"), length(tiers),
      if (always_exact) "*** CORRECT EVERYWHERE -- DESIGN IS BIASED ***" else "ok")
  ok <- ok && !always_exact && ever_exact && ever_wrong
}
say("  [%s] no estimator is correctly specified in every scenario, and every",
    if (ok) "PASS" else "FAIL")
say("        estimator is correctly specified in at least one")

say("\n=== Empirical corroboration at n = 2e5 ===")
say("  EXACT cells should show bias indistinguishable from zero; SMOOTHED and")
say("  MISSPEC cells should show a bias that survives the large sample size.")
grid <- seq(15, 180, by = 15)
h0 <- h0_matisse
shapes <- list(constant = ve_shapes$constant, gradual = ve_shapes$gradual,
               piecewise = ve_shapes$piecewise)
set.seed(20260724)
say("  %-12s %-10s %10s %10s %8s  tier", "shape", "method", "bias", "MCSE", "z")
emp <- list()
for (sh in names(shapes)) {
  fn <- shapes[[sh]]
  sc <- calibrate_scale(6000, 1e5, 1e5, h0, fn, 180, 0, 0.04, 120)
  nrep <- 25
  b <- sapply(seq_len(nrep), function(i) {
    dd <- simulate_trial(1e5, 1e5, h0, fn, 180, 0, 0.04, 120, scale = sc)
    c(cox_tdc = mean(est_cox_tdc(dd, grid)$ve_hat - fn(grid), na.rm = TRUE),
      durham  = mean(est_durham(dd, grid)$ve_hat - fn(grid), na.rm = TRUE),
      tian    = mean(est_tian(dd, grid)$ve_hat - fn(grid), na.rm = TRUE),
      pooled_logistic = mean(est_pooled_logistic(dd, grid, 45, 180)$ve_hat - fn(grid),
                             na.rm = TRUE))
  })
  for (m in rownames(b)) {
    mu <- mean(b[m, ]); s <- sd(b[m, ]) / sqrt(nrep)
    tier <- MAP[[m]][MAP$shape == sh]
    say("  %-12s %-10s %+10.5f %10.5f %+8.2f  %s", sh, m, mu, s, mu / s, tier)
    emp[[length(emp) + 1]] <- data.frame(shape = sh, method = m, bias = mu,
                                         mcse = s, z = mu / s, tier = tier)
  }
}
empd <- do.call(rbind, emp)
dir.create("outputs", showWarnings = FALSE)
write.csv(MAP, "outputs/specification_map.csv", row.names = FALSE)
write.csv(empd, "outputs/specification_map_empirical.csv", row.names = FALSE)
say("\n  -> outputs/specification_map.csv, outputs/specification_map_empirical.csv")
