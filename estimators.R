# ---------------------------------------------------------------------------
# Estimators of VE_h(t), all behind one interface.
#
# Every estimator returns a data.frame with columns t, ve_hat, lo, hi and the
# attributes: converged (logical), time_sec (numeric), note (character),
# diagnostics (list). Grid points where the method cannot produce an estimate
# return NA -- never a silently dropped row, and never a fabricated value.
#
# SIGN CONVENTION (brief sec 4): signed_error(t) = ve_hat(t) - VE_true(t),
# positive = overestimate. VE = 1 - exp(eta) is DECREASING in eta, so the UPPER
# eta bound maps to the LOWER VE bound. Every interval below goes through
# ve_from_eta(), which is the single place that inversion can be got wrong, and
# tests/testthat/test-inversion.R asserts it for all five methods.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(survival)
  library(splines)
})

# The one place eta -> VE happens. Note the deliberate swap of lo/hi.
ve_from_eta <- function(eta, se, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  eta_lo <- eta - z * se
  eta_hi <- eta + z * se
  list(ve_hat = 1 - exp(eta),
       lo     = 1 - exp(eta_hi),   # upper eta -> LOWER VE
       hi     = 1 - exp(eta_lo))
}

new_result <- function(t, ve_hat, lo, hi, converged, time_sec,
                       note = "", diagnostics = list()) {
  out <- data.frame(t = t, ve_hat = ve_hat, lo = lo, hi = hi)
  attr(out, "converged")   <- converged
  attr(out, "time_sec")    <- time_sec
  attr(out, "note")        <- note
  attr(out, "diagnostics") <- diagnostics
  out
}

fail_result <- function(grid, time_sec, note) {
  new_result(grid, NA_real_, NA_real_, NA_real_, FALSE, time_sec, note)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# 0. Fast risk-set core for a SINGLE BINARY covariate
# ---------------------------------------------------------------------------
# The Cox partial likelihood with one binary covariate needs, at each event time,
# only the number at risk in each arm. Writing exp(eta*z) with z in {0,1}:
#
#   S0(t) = n0(t) + n1(t)*exp(eta),   S1(t) = S2(t) = n1(t)*exp(eta)
#   zbar(t) = n1(t)e^eta / (n0(t) + n1(t)e^eta) = p(t)
#   S2/S0 - zbar^2 = p(1 - p)
#
# So every risk-set sum is a reverse cumulative count -- no stacked expansion,
# no per-event loop over subjects. Cox-TDC drops from ~7.4 s to milliseconds and
# Tian from ~72 s to well under a second, which is the difference between a
# feasible study and a 40-hour one. Correctness is not taken on trust:
# tests/testthat/test-fastcox.R asserts agreement with coxph(ties = "breslow")
# to 1e-6 on the coefficient and its standard error.
#
# Ties are handled by the Breslow approximation. Event times here are continuous,
# so exact ties are numerically negligible; the tie handling and a sensitivity
# check against Efron are reported rather than assumed away.

make_risk_core <- function(d) {
  o <- order(d$time)
  tt <- d$time[o]; zz <- d$arm[o]; st <- d$status[o]
  n <- length(tt)
  # number at risk (time >= t_j) in each arm, at each event time
  cum0 <- cumsum(zz == 0L); cum1 <- cumsum(zz == 1L)
  n0_tot <- cum0[n]; n1_tot <- cum1[n]
  ev <- which(st == 1L)
  # first index with time >= tt[j]  ->  at-risk = total - (count with time < t)
  lt <- findInterval(tt[ev] - 1e-12, tt)
  n0 <- n0_tot - c(0, cum0)[lt + 1L]
  n1 <- n1_tot - c(0, cum1)[lt + 1L]
  list(t_ev = tt[ev], z_ev = zz[ev], n0 = n0, n1 = n1, n_ev = length(ev))
}

# p(eta) at each event time, given the core
.core_p <- function(core, eta) {
  e <- exp(eta)
  a <- core$n1 * e
  a / (core$n0 + a)
}

# ---------------------------------------------------------------------------
# 1. Cox with a time-dependent covariate, g(t) = log t
# ---------------------------------------------------------------------------
# eta(t) = b1 + b2*log(t); delta method on the linear predictor.
# EXACT iff log(1 - VE_h(t)) is affine in log t; misspecified otherwise.

est_cox_tdc <- function(d, grid, maxit = 50, tol = 1e-10) {
  t0 <- proc.time()[["elapsed"]]
  core <- make_risk_core(d)
  if (core$n_ev < 3)
    return(fail_result(grid, proc.time()[["elapsed"]] - t0,
                       sprintf("only %d events", core$n_ev)))
  lg <- log(pmax(core$t_ev, 1e-8))
  X <- cbind(1, lg)
  b <- c(0, 0)
  conv <- FALSE
  for (it in seq_len(maxit)) {
    p <- .core_p(core, as.vector(X %*% b))
    r <- core$z_ev - p
    U <- crossprod(X, r)
    W <- p * (1 - p)
    I <- crossprod(X * W, X)
    step <- try(solve(I, U), silent = TRUE)
    if (inherits(step, "try-error") || any(!is.finite(step))) break
    b <- b + as.vector(step)
    if (any(!is.finite(b)) || max(abs(b)) > 50) break
    if (max(abs(step)) < tol) { conv <- TRUE; break }
  }
  el <- proc.time()[["elapsed"]] - t0
  if (!conv)
    return(fail_result(grid, el, "Cox-TDC Newton did not converge"))
  p <- .core_p(core, as.vector(X %*% b))
  I <- crossprod(X * (p * (1 - p)), X)
  V <- try(solve(I), silent = TRUE)
  if (inherits(V, "try-error"))
    return(fail_result(grid, el, "singular information matrix"))

  Xg  <- cbind(1, log(pmax(grid, 1e-8)))
  eta <- as.vector(Xg %*% b)
  se  <- sqrt(pmax(0, rowSums((Xg %*% V) * Xg)))
  v <- ve_from_eta(eta, se)
  # Wald test of no waning: H0 b2 = 0. This is the type I error / power test.
  z2 <- b[2] / sqrt(V[2, 2])
  new_result(grid, v$ve_hat, v$lo, v$hi, TRUE, el, "",
             list(beta = b, se_beta = sqrt(diag(V)),
                  wald_z_waning = z2,
                  wald_p_waning = 2 * pnorm(-abs(z2)),
                  n_events = core$n_ev))
}

# ---------------------------------------------------------------------------
# 2. Durham -- smoothed scaled Schoenfeld residuals
# ---------------------------------------------------------------------------
# This reproduces survival:::plot.cox.zph EXACTLY, which is the canonical
# implementation and NOT the GCV smoothing spline that methodology text often
# describes. Three details decide whether the bands are right:
#
#   (a) ns(x, df = 4, intercept = TRUE) -- the intercept is INSIDE the basis.
#       Adding a separate one gives 5 df and badly under-smooths.
#   (b) The variance scale is zph$var[1,1], the Grambsch-Therneau variance from
#       the Cox fit -- NOT an OLS residual variance from the spline regression.
#       Substituting the latter inflates bands into spurious blow-ups.
#   (c) transform = "identity". The cox.zph DEFAULT is "km", which puts the
#       spline on a 1 - KM(t) axis that is RANDOM ACROSS REPLICATES, so the
#       effective smoothing would differ replicate to replicate and the grid
#       would need remapping every time. Identity keeps the spline on calendar
#       time, which is the scale the study evaluates on.
#
# The band is +/- 2 SE, matching the canonical implementation, so its nominal
# level is 95.45%, not 95%. Stated, not silently corrected.

est_durham <- function(d, grid, df = 4) {
  t0 <- proc.time()[["elapsed"]]
  fit <- try(coxph(Surv(time, status) ~ arm, data = d), silent = TRUE)
  if (inherits(fit, "try-error"))
    return(fail_result(grid, proc.time()[["elapsed"]] - t0, "coxph error"))
  zph <- try(cox.zph(fit, transform = "identity"), silent = TRUE)
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(zph, "try-error"))
    return(fail_result(grid, el, "cox.zph error"))

  xx <- zph$x
  yy <- as.vector(zph$y)
  keep <- is.finite(xx) & is.finite(yy)
  if (sum(keep) <= df + 1)
    return(fail_result(grid, el, sprintf("only %d usable residuals for df = %d",
                                         sum(keep), df)))

  xrange <- range(xx[keep])
  bas <- try(ns(xx[keep], df = df, intercept = TRUE), silent = TRUE)
  if (inherits(bas, "try-error"))
    return(fail_result(grid, el, "ns() basis failure"))

  qmat <- qr(bas)
  if (qmat$rank < df)
    return(fail_result(grid, el, sprintf("spline basis rank %d < df %d", qmat$rank, df)))

  pmat <- predict(bas, newx = grid)
  eta  <- as.vector(pmat %*% qr.coef(qmat, yy[keep]))
  bk    <- backsolve(qr.R(qmat)[1:df, 1:df], diag(df))
  seval <- rowSums((pmat %*% (bk %*% t(bk))) * pmat)
  se    <- sqrt(pmax(0, zph$var[1, 1] * seval))

  # +/- 2 SE to match the canonical implementation (nominal 95.45%)
  v <- list(ve_hat = 1 - exp(eta),
            lo     = 1 - exp(eta + 2 * se),
            hi     = 1 - exp(eta - 2 * se))
  outside <- grid < xrange[1] | grid > xrange[2]
  new_result(grid, v$ve_hat, v$lo, v$hi, TRUE, el,
             if (any(outside)) sprintf("%d grid points extrapolated beyond event support",
                                       sum(outside)) else "",
             list(realised_df = qmat$rank, gt_var = zph$var[1, 1],
                  zph_p = zph$table[1, "p"], extrapolated = outside,
                  n_resid = sum(keep)))
}

# ---------------------------------------------------------------------------
# 3. Tian -- kernel-weighted local partial likelihood
# ---------------------------------------------------------------------------
# Epanechnikov kernel on event times; scalar Newton solve at each grid point;
# robust (sandwich) local standard errors with squared kernel weights.
# Bandwidth by K-fold cross-validated partial likelihood, per Tian et al. (2005).

.epan <- function(u) ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)

# Local fit at one grid point. Fully vectorised via the binary-covariate core.
.tian_at <- function(t0, h, core, maxit = 40, tol = 1e-9) {
  w <- .epan((core$t_ev - t0) / h) / h
  use <- w > 0
  neff <- sum(use)
  if (neff < 3) return(list(beta = NA_real_, se = NA_real_, neff = neff))
  ww <- w[use]
  sub <- list(n0 = core$n0[use], n1 = core$n1[use])
  z <- core$z_ev[use]
  b <- 0
  conv <- FALSE
  for (it in seq_len(maxit)) {
    p <- .core_p(sub, b)
    U <- sum(ww * (z - p))
    I <- sum(ww * p * (1 - p))
    if (!is.finite(I) || I <= 1e-14) return(list(beta = NA_real_, se = NA_real_, neff = neff))
    step <- U / I
    b <- b + step
    if (!is.finite(b) || abs(b) > 20) return(list(beta = NA_real_, se = NA_real_, neff = neff))
    if (abs(step) < tol) { conv <- TRUE; break }
  }
  if (!conv) return(list(beta = NA_real_, se = NA_real_, neff = neff))
  p <- .core_p(sub, b)
  I <- sum(ww * p * (1 - p))
  J <- sum(ww^2 * (z - p)^2)          # robust local sandwich meat
  se <- if (I > 1e-14 && J >= 0) sqrt(J) / I else NA_real_
  list(beta = b, se = se, neff = neff)
}

# Cross-validated partial likelihood for bandwidth selection (Tian et al. 2005:
# minus log partial likelihood as prediction error).
.tian_cv <- function(d, h_grid, K = 5, grid_cv) {
  n <- nrow(d)
  fold <- sample(rep_len(seq_len(K), n))
  vapply(h_grid, function(h) {
    tot <- 0
    for (k in seq_len(K)) {
      tr_core <- make_risk_core(d[fold != k, ])
      if (tr_core$n_ev < 5) return(-Inf)
      bs <- vapply(grid_cv, function(t0) .tian_at(t0, h, tr_core)$beta, 0)
      if (all(is.na(bs))) return(-Inf)
      bs <- approx(grid_cv[!is.na(bs)], bs[!is.na(bs)], xout = grid_cv, rule = 2)$y
      te_core <- make_risk_core(d[fold == k, ])
      if (te_core$n_ev == 0) next
      bt <- approx(grid_cv, bs, xout = te_core$t_ev, rule = 2)$y
      ok <- is.finite(bt)
      if (!any(ok)) next
      e <- exp(bt[ok])
      tot <- tot + sum(bt[ok] * te_core$z_ev[ok] -
                         log(te_core$n0[ok] + te_core$n1[ok] * e))
    }
    tot
  }, 0)
}

est_tian <- function(d, grid, h_grid = NULL, cv_folds = 5) {
  t0 <- proc.time()[["elapsed"]]
  core <- make_risk_core(d)
  if (core$n_ev < 10)
    return(fail_result(grid, proc.time()[["elapsed"]] - t0,
                       sprintf("only %d events", core$n_ev)))
  span <- diff(range(grid))
  # The search must extend past h = span, or the boundary itself constrains the
  # answer. The Gate 4 pilot selected the top of an earlier grid (0.75 * span)
  # in 62-76% of replicates, which means cross-validation wanted MORE smoothing
  # than the grid allowed and the reported bandwidths were an artefact of the
  # truncation. At h >= ~2 * span the Epanechnikov window covers the whole
  # follow-up and the estimator degenerates to a global constant fit -- that is
  # a legitimate answer for CV to give at these event counts, and the study
  # should be able to observe it rather than being prevented from doing so.
  if (is.null(h_grid)) h_grid <- span * c(0.15, 0.25, 0.40, 0.60, 0.90, 1.30, 2.00)
  grid_cv <- seq(min(grid), max(grid), length.out = 8)
  cv <- try(.tian_cv(d, h_grid, cv_folds, grid_cv), silent = TRUE)
  h <- if (inherits(cv, "try-error") || all(!is.finite(cv))) {
    h_grid[ceiling(length(h_grid) / 2)]
  } else h_grid[which.max(cv)]
  at_boundary <- isTRUE(all.equal(h, h_grid[1])) ||
    isTRUE(all.equal(h, h_grid[length(h_grid)]))

  fits <- lapply(grid, function(g) .tian_at(g, h, core))
  eta <- vapply(fits, function(f) f$beta, 0)
  se  <- vapply(fits, function(f) f$se, 0)
  el <- proc.time()[["elapsed"]] - t0
  ok <- is.finite(eta) & is.finite(se)
  v <- ve_from_eta(eta, se)
  v$ve_hat[!ok] <- NA; v$lo[!ok] <- NA; v$hi[!ok] <- NA
  new_result(grid, v$ve_hat, v$lo, v$hi, any(ok), el,
             if (all(ok)) "" else sprintf("%d/%d grid points not estimable",
                                          sum(!ok), length(grid)),
             list(bandwidth = h, bandwidth_at_boundary = at_boundary,
                  cv = cv, n_eff = vapply(fits, function(f) f$neff, 0)))
}

# ---------------------------------------------------------------------------
# 4. Pooled logistic with time-varying effects
# ---------------------------------------------------------------------------
# Discrete-time grouped survival on an interval grid, with interval-specific
# treatment effects.
#
# COMPLEMENTARY LOG-LOG link, not logit. Under cloglog the coefficient IS the
# log hazard ratio for grouped survival exactly, so no rare-event approximation
# is needed anywhere and eta feeds straight into ve_from_eta(). With logit one
# would be estimating a log odds ratio and relying on events being rare.
#
# The interval grid is deliberately COARSER than the DGP's 30-day hazard blocks
# so the model is not saturated against the truth (see brief sec 3.2).

est_pooled_logistic <- function(d, grid, interval_width = 45, horizon = NULL) {
  t0 <- proc.time()[["elapsed"]]
  if (is.null(horizon)) horizon <- max(grid)
  brk <- seq(0, horizon, by = interval_width)
  if (tail(brk, 1) < horizon) brk <- c(brk, horizon)
  ni <- length(brk) - 1L

  # person-interval expansion
  rows <- lapply(seq_len(ni), function(k) {
    lo <- brk[k]; hi <- brk[k + 1]
    at_risk <- d$time > lo
    if (!any(at_risk)) return(NULL)
    sub <- d[at_risk, ]
    ev <- as.integer(sub$status == 1L & sub$time <= hi)
    expo <- pmin(sub$time, hi) - lo
    data.frame(interval = factor(k, levels = seq_len(ni)), arm = sub$arm,
               event = ev, logexp = log(pmax(expo, 1e-6)))
  })
  pl <- do.call(rbind, rows[!vapply(rows, is.null, TRUE)])
  el0 <- proc.time()[["elapsed"]] - t0
  if (is.null(pl) || nrow(pl) == 0 || sum(pl$event) < 2)
    return(fail_result(grid, el0, "no usable person-intervals"))

  fit <- try(suppressWarnings(
    glm(event ~ interval + interval:arm - 1 + offset(logexp),
        family = binomial(link = "cloglog"), data = pl)), silent = TRUE)
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(fit, "try-error"))
    return(fail_result(grid, el, "glm error"))
  if (!fit$converged)
    return(fail_result(grid, el, "glm did not converge"))

  cf <- coef(fit); V <- vcov(fit)
  arm_idx <- grep(":arm$", names(cf))
  if (length(arm_idx) != ni)
    return(fail_result(grid, el, "interval effects not fully identified"))
  # separation: an interval with zero events in one arm gives |coef| huge
  sep <- !is.finite(cf[arm_idx]) | abs(cf[arm_idx]) > 15

  gi <- pmin(ni, pmax(1L, findInterval(grid, brk, rightmost.closed = TRUE)))
  eta <- cf[arm_idx][gi]
  se  <- sqrt(pmax(0, diag(V)[arm_idx][gi]))
  bad <- sep[gi] | !is.finite(eta) | !is.finite(se)
  v <- ve_from_eta(eta, se)
  v$ve_hat[bad] <- NA; v$lo[bad] <- NA; v$hi[bad] <- NA
  new_result(grid, v$ve_hat, v$lo, v$hi, TRUE, el,
             if (any(sep)) sprintf("%d interval(s) separated (zero events in an arm)",
                                   sum(sep)) else "",
             list(n_intervals = ni, interval_width = interval_width,
                  separated = sep, breaks = brk))
}
