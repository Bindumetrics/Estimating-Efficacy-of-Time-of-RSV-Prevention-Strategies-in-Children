# ---------------------------------------------------------------------------
# GATE 1 -- DGP correctness. Every check reports a NUMBER.
# Run: Rscript scripts/gate1_dgp.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(survival) })
source("R/dgp.R")

set.seed(20260721)
PASS <- TRUE
res <- list()
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
check <- function(label, value, tol, unit = "") {
  ok <- is.finite(value) && abs(value) <= tol
  PASS <<- PASS && ok
  res[[label]] <<- value
  say("  [%s] %-46s %12.6g %s (tol %.4g)", if (ok) "PASS" else "FAIL",
      label, value, unit, tol)
  invisible(ok)
}

h0 <- h0_matisse
horizon <- 180

# ---------------------------------------------------------------------------
say("\n=== 1.1 Analytic survival check (control arm, homogeneous) ===")
# Large control arm, no censoring at all; compare empirical KM to exp(-H(t)).
n_big <- 400000
d <- simulate_trial(n_c = n_big, n_v = 0, h0 = h0, ve_id_fn = ve_shapes$constant,
                    horizon = 1e9, theta = 0, ltfu_prob = 0, accrual_days = 0)
km <- survfit(Surv(time, status) ~ 1, data = d)
tgrid <- seq(15, 360, by = 15)
s_emp <- summary(km, times = tgrid, extend = TRUE)$surv
s_ana <- exp(-pc_integral(h0, tgrid))
say("  n = %d, events = %d", n_big, sum(d$status))
say("  max |KM - analytic| over t in [15,360] by 15:")
check("survival.max_abs_dev", max(abs(s_emp - s_ana)), tol = 0.004)
# Monte Carlo yardstick: pointwise SE of KM at the worst point
se_km <- summary(km, times = tgrid, extend = TRUE)$std.err
say("  (max pointwise KM standard error = %.6f, so the deviation above should",
    max(se_km))
say("   be of that order -- it is %.2f x the max SE)",
    max(abs(s_emp - s_ana)) / max(se_km))

# ---------------------------------------------------------------------------
say("\n=== 1.2 Hazard-ratio recovery (sampler validation) ===")
# NOTE: under theta = 0 this is near-tautological by construction (see sec 3.0
# of the brief). It validates the sampler and the inversion, not the design.
#
# Scored as a STANDARDISED deviation, not an absolute one. An absolute tolerance
# on a ratio is the wrong statistic: Var(log HR_hat) = 1/D0 + 1/D1, so the
# sampling error of HR_hat scales WITH HR and inversely with sqrt(events). For
# the `rapid` shape at late t (HR ~ 0.95, ~1200 events/bin) the analytic SE is
# ~0.039, so an absolute deviation of 0.08 is ordinary noise, while for the
# `constant` shape (HR = 0.30) the same 0.08 would be a real defect. z is
# comparable across bins and shapes; raw deviations are reported alongside.
zmax_all <- 0
nbin_tot <- 0
for (shape in c("constant", "gradual", "rapid", "piecewise")) {
  fn <- ve_shapes[[shape]]
  n_arm <- 300000
  dd <- simulate_trial(n_arm, n_arm, h0, fn, horizon = 1e9, theta = 0,
                       ltfu_prob = 0, accrual_days = 0)
  brk <- seq(0, 360, by = 20)
  mid <- head(brk, -1) + diff(brk) / 2
  st <- t(sapply(seq_along(mid), function(k) {
    lo <- brk[k]; hi <- brk[k + 1]
    inb <- dd$status == 1L & dd$time > lo & dd$time <= hi
    d0 <- sum(inb & dd$arm == 0L); d1 <- sum(inb & dd$arm == 1L)
    py <- tapply(pmin(dd$time, hi) - pmin(dd$time, lo), dd$arm, sum)
    c(d0 = d0, d1 = d1, hr = (d1 / py[["1"]]) / (d0 / py[["0"]]))
  }))
  hr_true <- 1 - fn(mid)
  z <- (log(st[, "hr"]) - log(hr_true)) / sqrt(1 / st[, "d0"] + 1 / st[, "d1"])
  dev <- abs(st[, "hr"] - hr_true)
  say("  %-10s  max|z| = %5.2f   mean z = %+5.2f   sd z = %4.2f   (max abs dev %.4f, median events/bin %d)",
      shape, max(abs(z)), mean(z), sd(z), max(dev), as.integer(median(st[, "d1"])))
  res[[paste0("hr.", shape, ".max_abs_z")]] <- max(abs(z))
  res[[paste0("hr.", shape, ".max_abs_dev")]] <- max(dev)
  zmax_all <- max(zmax_all, max(abs(z))); nbin_tot <- nbin_tot + length(z)
}
# Bonferroni-style yardstick: max |z| over N independent standard normals has
# expected max ~ qnorm(1 - 0.05/(2N)). Exceeding that is evidence of real bias.
crit <- qnorm(1 - 0.05 / (2 * nbin_tot))
ok12 <- zmax_all <= crit
PASS <- PASS && ok12
say("  [%s] max|z| = %.2f over %d bins vs critical value %.2f (5%% family-wise)",
    if (ok12) "PASS" else "FAIL", zmax_all, nbin_tot, crit)

# ---------------------------------------------------------------------------
say("\n=== 1.3 Inversion accuracy: numeric vs closed form ===")
# The piecewise VE shape admits a closed-form treated-arm inverse. Drive both
# paths from the SAME uniform draws and compare event times directly.
brks <- c(100, 135, 170); lvls <- c(0.80, 0.65, 0.50, 0.35)
fn_pw <- ve_piecewise(brks, lvls)
set.seed(99); e <- rexp(200000)
# closed form
cuts <- sort(unique(c(h0$cuts, brks)))
vals <- pc_eval(h0, cuts) * (1 - lvls[findInterval(cuts, c(0, brks))])
t_exact <- pc_integral_inverse(list(cuts = cuts, vals = vals), e)
# numeric path (the one used for smooth shapes)
grid <- seq(0, 800, by = 0.25)
hv <- pc_eval(h0, grid) * (1 - fn_pw(grid))
cum <- c(0, cumsum((head(hv, -1) + tail(hv, -1)) / 2 * 0.25))
t_num <- approx(cum, grid, xout = e, rule = 1)$y
keep <- is.finite(t_exact) & is.finite(t_num) & t_exact <= 400
say("  compared on %d draws with T <= 400 d", sum(keep))
check("inversion.max_abs_days", max(abs(t_exact[keep] - t_num[keep])), tol = 0.15, "days")
check("inversion.median_abs_days", median(abs(t_exact[keep] - t_num[keep])), tol = 0.15, "days")

# ---------------------------------------------------------------------------
say("\n=== 1.4 Event-count calibration ===")
targets <- c(93, 200, 304)
horizons <- c(180, 180, 360)
ok14 <- TRUE
for (i in seq_along(targets)) {
  tg <- targets[i]; hz <- horizons[i]
  sc <- calibrate_scale(tg, 3585, 3563, h0, ve_shapes$gradual, hz,
                        theta = 0, ltfu_prob = 0.04, accrual_days = 120)
  nrep <- 400
  ev <- replicate(nrep, sum(simulate_trial(3585, 3563, h0, ve_shapes$gradual, hz,
                                           theta = 0, ltfu_prob = 0.04,
                                           accrual_days = 120, scale = sc)$status))
  pe <- 100 * (mean(ev) - tg) / tg
  mcse <- 100 * (sd(ev) / sqrt(nrep)) / tg
  say("  target %3d @ %3d d: scale = %.4f, mean observed = %.1f (MCSE %.1f), error = %+.2f%% (MCSE %.2f%%)",
      tg, hz, sc, mean(ev), sd(ev) / sqrt(nrep), pe, mcse)
  res[[paste0("calib.", tg)]] <- pe
  ok14 <- ok14 && abs(pe) < 2
}
PASS <- PASS && ok14
say("  [%s] all realised event counts within 2%% of target", if (ok14) "PASS" else "FAIL")

# ---------------------------------------------------------------------------
say("\n=== 1.5 Censoring check ===")
sc <- calibrate_scale(200, 3585, 3563, h0, ve_shapes$gradual, 180, 0, 0.04, 0)
nrep <- 300
frac <- replicate(nrep, {
  dd <- simulate_trial(3585, 3563, h0, ve_shapes$gradual, 180, theta = 0,
                       ltfu_prob = 0.04, accrual_days = 0, scale = sc)
  mean(is.finite(dd$t_cens) & dd$t_cens < 180)
})
say("  nominal LTFU by 180 d = 4.00%%, realised = %.2f%% (MCSE %.3f%%)",
    100 * mean(frac), 100 * sd(frac) / sqrt(nrep))
check("censoring.ltfu_pct_error", 100 * mean(frac) - 4, tol = 0.25, "pp")

# independence: event-time distribution should not differ between subjects
# censored early vs late. Compare LATENT event times across censoring halves.
dd <- simulate_trial(60000, 60000, h0, ve_shapes$gradual, 180, theta = 0,
                     ltfu_prob = 0.20, accrual_days = 0, scale = sc)
cens <- dd[is.finite(dd$t_cens) & dd$t_cens < 180, ]
half <- cens$t_cens < median(cens$t_cens)
kt <- ks.test(cens$t_event[half], cens$t_event[!half])
say("  independence of censoring: KS on latent event times, early- vs late-censored")
say("  D = %.5f, p = %.3f (n = %d vs %d) -- expect p not small",
    kt$statistic, kt$p.value, sum(half), sum(!half))
res[["censoring.ks_p"]] <- kt$p.value
PASS <- PASS && kt$p.value > 0.01
say("  [%s] no evidence of dependence at the 1%% level", if (kt$p.value > 0.01) "PASS" else "FAIL")

# ---------------------------------------------------------------------------
say("\n=== 1.6 Degenerate-case check: exit-time pile-up ===")
say("  Fraction of ALL subjects exiting exactly at the horizon:")
for (cfg in list(list(l = 0, a = 0,   lab = "admin-only, no accrual (degenerate)"),
                 list(l = 0, a = 120, lab = "staggered entry, no LTFU"),
                 list(l = 0.04, a = 120, lab = "reference: 4% LTFU + accrual"),
                 list(l = 0.20, a = 0,   lab = "20% LTFU, no accrual"))) {
  dd <- simulate_trial(3585, 3563, h0, ve_shapes$gradual, 180, theta = 0,
                       ltfu_prob = cfg$l, accrual_days = cfg$a, scale = sc)
  pile <- mean(dd$time >= 180 - 1e-9)
  n_int <- length(unique(round(dd$time[dd$time < 180 - 1e-9], 1)))
  say("    %-38s pile-up = %5.1f%%   distinct interior exit times = %d",
      cfg$lab, 100 * pile, n_int)
  res[[paste0("pileup.", cfg$lab)]] <- pile
}
say("  -> the first row is the configuration that degenerates Tian's")
say("     observation-time kernel weights (brief sec 9). Reported, not hidden.")

# ---------------------------------------------------------------------------
say("\n=== 1.7 Frailty: VE_h(t) must differ from VE_id(t) ===")
# Constant individual VE, theta > 0. The MARGINAL hazard ratio must wane.
tg <- seq(15, 360, by = 15)
for (theta in c(0, 0.5, 1.0)) {
  sc2 <- calibrate_scale(304, 3585, 3563, h0, ve_shapes$constant, 360, theta, 0.04, 0)
  h0s <- list(cuts = h0$cuts, vals = h0$vals * sc2)
  vh <- ve_h_true(tg, h0s, ve_shapes$constant, theta)
  say("  theta = %.1f: VE_id = %.3f constant; VE_h at 30/180/360 d = %.4f / %.4f / %.4f  (drift %.4f)",
      theta, 0.70, vh[2], vh[12], vh[24], vh[2] - vh[24])
  res[[paste0("frailty.drift.theta", theta)]] <- vh[2] - vh[24]
}
say("  -> nonzero drift at theta > 0 is depletion of susceptibles: apparent")
say("     waning with NO individual waning. No estimator in the panel targets VE_id.")

# ---------------------------------------------------------------------------
say("\n=== 1.8 External consistency: simulated vs published cumulative VE ===")
# THE check with teeth. Published MATISSE final-analysis cumulative VE, any-severity.
pub_t  <- c(90, 120, 150, 180, 210, 240, 270, 360)
pub_ve <- c(57.6, 54.5, 50.0, 49.2, 43.8, 39.7, 35.0, 33.0) / 100
pub_lo <- c(31.3, 33.2, 30.3, 31.4, 25.6, 21.3, 16.1, 15.2) / 100
pub_hi <- c(74.6, 69.5, 64.5, 62.8, 57.7, 54.1, 49.9, 47.1) / 100
# Sampling SE of each published VE, recovered from its reported 95% CI. The
# targets are themselves estimates from ~300 events, with SEs of 8-11pp; scoring
# the simulation against them to a flat 2pp would be demanding agreement far
# tighter than the data can support. Deviations are therefore reported BOTH in
# percentage points and in units of the published SE.
pub_se <- (pub_hi - pub_lo) / (2 * qnorm(0.975))

# Fit the log-quadratic family (declared a priori in R/dgp.R; one degree above
# the Cox-TDC form, so Cox-TDC stays misspecified against it by construction).
obj <- function(p) {
  fn <- ve_logquad(p[1], p[2], p[3])
  sc <- calibrate_scale(304, 3585, 3563, h0, fn, 360, 0, 0.04, 0)
  h0s <- list(cuts = h0$cuts, vals = h0$vals * sc)
  F_C <- 1 - exp(-pc_integral(h0s, pub_t))
  F_V <- 1 - exp(-cumhaz_treated(pub_t, h0s, fn))
  sum(((1 - F_V / F_C - pub_ve) / pub_se)^2)   # weighted by published precision
}
op <- optim(c(log(0.4), 0.2, 0.05), obj, method = "Nelder-Mead",
            control = list(maxit = 2000, reltol = 1e-12))
fit <- ve_logquad(op$par[1], op$par[2], op$par[3])
sc <- calibrate_scale(304, 3585, 3563, h0, fit, 360, 0, 0.04, 0)
h0s <- list(cuts = h0$cuts, vals = h0$vals * sc)
F_C <- 1 - exp(-pc_integral(h0s, pub_t))
F_V <- 1 - exp(-cumhaz_treated(pub_t, h0s, fit))
sim_ve <- 1 - F_V / F_C
say("  fitted VE_id (log-quadratic): a = %.4f, b = %.4f, c = %.4f",
    op$par[1], op$par[2], op$par[3])
say("   day  published (95%% CI)        SE     simulated    diff      z")
for (k in seq_along(pub_t))
  say("  %4d    %5.1f%% (%4.1f-%4.1f)   %4.1fpp    %5.1f%%   %+5.2fpp  %+5.2f",
      pub_t[k], 100 * pub_ve[k], 100 * pub_lo[k], 100 * pub_hi[k],
      100 * pub_se[k], 100 * sim_ve[k], 100 * (sim_ve[k] - pub_ve[k]),
      (sim_ve[k] - pub_ve[k]) / pub_se[k])
zmax <- max(abs((sim_ve - pub_ve) / pub_se))
say("  max |deviation| = %.2f pp = %.2f published SE",
    100 * max(abs(sim_ve - pub_ve)), zmax)
res[["published.max_abs_dev_pp"]] <- 100 * max(abs(sim_ve - pub_ve))
ok18 <- zmax <= 1.0
PASS <- PASS && ok18
res[["published.max_abs_z"]] <- zmax
say("  [%s] every simulated cumulative VE within 1 published SE of its target",
    if (ok18) "PASS" else "FAIL")
say("  instantaneous VE_id at 90/180/360 d = %.3f / %.3f / %.3f",
    fit(90), fit(180), fit(360))
saveRDS(list(par = op$par, scale = sc), "outputs/environment/applied_anchor_fit.rds")
say("  -> note the instantaneous curve sits BELOW cumulative VE at every t, as")
say("     it must: cumulative VE is a weighted average that lags the decline.")

# ---------------------------------------------------------------------------
say("\n=== GATE 1 SUMMARY ===")
say("  Overall: %s", if (PASS) "PASS" else "FAIL")
dir.create("outputs/environment", showWarnings = FALSE, recursive = TRUE)
saveRDS(res, "outputs/environment/gate1_results.rds")
write.csv(data.frame(check = names(res), value = unlist(res)),
          "outputs/environment/gate1_results.csv", row.names = FALSE)
say("  Numbers written to outputs/environment/gate1_results.csv")
