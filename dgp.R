# ---------------------------------------------------------------------------
# Data-generating mechanism for the VE(t) simulation study.
#
# Individual-level survival data. Piecewise-constant baseline hazard on 30-day
# blocks; multiplicative gamma frailty; two independent censoring mechanisms.
#
# Estimand (SIMULATION_PROMPT.md sec 0.7, 3.0): the POPULATION hazard-ratio VE,
#   VE_h(t) = 1 - lambda_V(t)/lambda_C(t)
# where lambda are MARGINAL arm hazards. Under frailty this is not equal to the
# individual-level VE_id(t) that the DGP is parameterised by.
# ---------------------------------------------------------------------------

# --- piecewise-constant hazard utilities -----------------------------------

# A piecewise-constant function on [0, Inf) is stored as list(cuts, vals) where
# cuts = c(0, t1, t2, ..., tK) are the left endpoints and vals[k] applies on
# [cuts[k], cuts[k+1]).

pc_eval <- function(pc, t) {
  idx <- findInterval(t, pc$cuts, rightmost.closed = FALSE)
  idx[idx < 1L] <- 1L
  pc$vals[idx]
}

# Integral from 0 to t of a piecewise-constant function. Vectorised over t.
pc_integral <- function(pc, t) {
  cuts <- pc$cuts
  vals <- pc$vals
  cum_at_cut <- c(0, cumsum(vals[-length(vals)] * diff(cuts)))
  idx <- findInterval(t, cuts, rightmost.closed = FALSE)
  idx[idx < 1L] <- 1L
  cum_at_cut[idx] + vals[idx] * (t - cuts[idx])
}

# Inverse of pc_integral: smallest t with integral == target. Closed form,
# because a piecewise-constant hazard has a piecewise-LINEAR cumulative hazard.
# Returns Inf when target exceeds the total integral over the final block, which
# for an open-ended final block cannot happen unless vals[K] == 0.
pc_integral_inverse <- function(pc, target) {
  cuts <- pc$cuts
  vals <- pc$vals
  cum_at_cut <- c(0, cumsum(vals[-length(vals)] * diff(cuts)))
  # which block does the target land in?
  idx <- findInterval(target, cum_at_cut, rightmost.closed = FALSE)
  idx[idx < 1L] <- 1L
  out <- cuts[idx] + (target - cum_at_cut[idx]) / vals[idx]
  out[vals[idx] <= 0] <- Inf
  out
}

# --- true VE_id(t) shapes ---------------------------------------------------
#
# Each returns a function of t giving the INDIVIDUAL-level hazard-scale VE.
# None of these is the fitted family of any single competitor; the piecewise
# shape's breakpoints are deliberately off the 30-day hazard blocks.

ve_constant <- function(ve0 = 0.70) {
  force(ve0)
  function(t) rep(ve0, length(t))
}

# Slow monotone decay toward a plateau. Not an exponential decay: uses a
# logistic-in-log-time form, deliberately outside the Bayesian candidate set.
ve_gradual <- function(ve0 = 0.80, ve_inf = 0.45, t_mid = 200, slope = 2.2) {
  function(t) {
    w <- 1 / (1 + (pmax(t, 1e-8) / t_mid)^slope)
    ve_inf + (ve0 - ve_inf) * w
  }
}

# Steep early decay flattening to near zero. Same family, harsher tuning.
ve_rapid <- function(ve0 = 0.85, ve_inf = 0.05, t_mid = 75, slope = 3.0) {
  function(t) {
    w <- 1 / (1 + (pmax(t, 1e-8) / t_mid)^slope)
    ve_inf + (ve0 - ve_inf) * w
  }
}

# Step changes at 100/135/170 d -- deliberately NOT on the 30-day hazard blocks,
# so no interval estimator's grid coincides with the truth's breakpoints.
ve_piecewise <- function(breaks = c(100, 135, 170),
                         levels = c(0.80, 0.65, 0.50, 0.35)) {
  stopifnot(length(levels) == length(breaks) + 1L)
  pc <- list(cuts = c(0, breaks), vals = levels)
  function(t) pc_eval(pc, t)
}

# Applied-anchor family: log(1 - VE_id(t)) quadratic in log(1 + t/30).
#
# Chosen a priori and flexible enough to track a real trial's cumulative-VE
# trajectory, which the logistic families above are not.
#
# NOTE on its relation to Cox-TDC, which fits log(1 - VE(t)) = b1 + b2*log(t).
# The basis here is log(1 + t/30), and log(1 + t/30) is NOT an affine function of
# log(t) -- the shift by 1 breaks it. So Cox-TDC is misspecified against this
# family for EVERY (a, b, c), including c = 0, not merely when c != 0. That is
# stronger than the design needs, but it means ve_logquad(c = 0) must NOT be used
# as the "Cox-TDC is exactly specified" truth; use ve_logt() below for that.
ve_logquad <- function(a, b, c) {
  function(t) {
    x <- log1p(pmax(t, 0) / 30)
    1 - exp(a + b * x + c * x^2)
  }
}

# The Cox-TDC model's OWN family: log(1 - VE(t)) affine in log t. Used only to
# construct the truth under which Cox-TDC is exactly specified, for the interval
# calibration check in Gate 2.4. Not part of the scenario grid -- putting an
# estimator's own family into the grid as a headline scenario is the circularity
# the brief forbids.
ve_logt <- function(a, b) {
  function(t) 1 - exp(a + b * log(pmax(t, 1e-8)))
}

ve_shapes <- list(
  constant  = ve_constant(),
  gradual   = ve_gradual(),
  rapid     = ve_rapid(),
  piecewise = ve_piecewise()
)

# --- marginal (population) hazard-ratio VE under gamma frailty ---------------
#
# With U ~ Gamma(shape = 1/theta, scale = theta), mean 1 and variance theta,
# the marginal hazard is lambda_a(t) = h_a(t) / (1 + theta * H_a(t)), so
#
#   1 - VE_h(t) = {1 - VE_id(t)} * (1 + theta*H_C(t)) / (1 + theta*H_V(t))
#
# theta = 0 gives VE_h == VE_id exactly.

ve_h_true <- function(t, h0, ve_id_fn, theta) {
  ve_id <- ve_id_fn(t)
  if (theta <= 0) return(ve_id)
  H_C <- pc_integral(h0, t)
  H_V <- cumhaz_treated(t, h0, ve_id_fn)
  1 - (1 - ve_id) * (1 + theta * H_C) / (1 + theta * H_V)
}

# Cumulative treated hazard H_V(t) = int_0^t h0(s) * {1 - VE_id(s)} ds.
# Exact when VE_id is piecewise-constant on the same cuts as h0; otherwise
# integrated numerically on a fine grid (trapezoid on 0.25-day steps, which at
# these hazard magnitudes is accurate to ~1e-10 relative).
#
# Integration uses the MIDPOINT rule, not the trapezoid. The integrand is
# piecewise-constant in its h0 factor, and all h0/VE knots are chosen to land on
# grid points, so no cell ever contains an interior discontinuity and midpoint is
# exact for that factor. The trapezoid rule straddles knots and mis-integrates
# each by (delta_h * step / 2), which at these hazard magnitudes showed up as
# ~0.1 day of inversion error -- small, but avoidable for free.
cumhaz_treated <- function(t, h0, ve_id_fn, step = 0.25) {
  tmax <- max(t, na.rm = TRUE)
  grid <- seq(0, tmax + step, by = step)
  mids <- head(grid, -1) + step / 2
  hv <- pc_eval(h0, mids) * (1 - ve_id_fn(mids))
  cum <- c(0, cumsum(hv * step))
  approx(grid, cum, xout = t, rule = 2)$y
}

# Closed-form version, valid only when VE_id is piecewise-constant on cuts that
# are a subset of a common refinement with h0. Used as the analytic reference in
# the Gate 1 inversion check.
cumhaz_treated_exact <- function(t, h0, ve_breaks, ve_levels) {
  cuts <- sort(unique(c(h0$cuts, ve_breaks)))
  vals <- pc_eval(h0, cuts) * (1 - ve_levels[findInterval(cuts, c(0, ve_breaks))])
  pc_integral(list(cuts = cuts, vals = vals), t)
}

# --- event-time sampling ----------------------------------------------------
#
# Inverse transform: draw E ~ Exp(1), solve H(T) = E/U for T.
# Control arm: closed form (piecewise-linear cumulative hazard).
# Treated arm: closed form when VE_id is piecewise-constant on the hazard blocks,
# otherwise monotone root-finding on a fine grid with linear interpolation.

sample_event_times <- function(n, h0, ve_id_fn, frailty, scale = 1,
                               exact_pc = NULL, step = 0.25, tmax = 800) {
  e <- rexp(n) / frailty
  h0s <- list(cuts = h0$cuts, vals = h0$vals * scale)
  if (is.null(ve_id_fn)) {
    return(pc_integral_inverse(h0s, e))
  }
  if (!is.null(exact_pc)) {
    cuts <- sort(unique(c(h0s$cuts, exact_pc$breaks)))
    vals <- pc_eval(h0s, cuts) *
      (1 - exact_pc$levels[findInterval(cuts, c(0, exact_pc$breaks))])
    return(pc_integral_inverse(list(cuts = cuts, vals = vals), e))
  }
  grid <- seq(0, tmax, by = step)
  mids <- head(grid, -1) + step / 2
  hv <- pc_eval(h0s, mids) * (1 - ve_id_fn(mids))
  cum <- c(0, cumsum(hv * step))
  # invert by interpolation; monotone because hv >= 0
  out <- approx(cum, grid, xout = e, rule = 1)$y
  out[is.na(out)] <- Inf   # E beyond total hazard over [0, tmax] -> no event
  out
}

# --- the generator ----------------------------------------------------------

#' Simulate one trial dataset.
#'
#' @param n_c,n_v arm sizes (control, vaccine)
#' @param h0 piecewise-constant baseline hazard, list(cuts, vals)
#' @param ve_id_fn function(t) giving individual-level hazard-scale VE
#' @param horizon administrative follow-up horizon in days
#' @param theta gamma frailty variance; 0 = homogeneous
#' @param ltfu_prob probability of loss to follow-up by `horizon` (exponential)
#' @param accrual_days width of the uniform staggered-entry window; 0 = off.
#' @param admin_complete_frac fraction of subjects who reach the full horizon
#'   before the data cutoff. MATISSE states that "data from 85% of the scheduled
#'   follow-up through 180 days were available", so 0.85 is the verified value.
#'
#'   Entry ~ U(0, A) and the cutoff sits at calendar time C = horizon + f*A, so
#'   available follow-up is min(horizon, horizon + f*A - entry). A subject with
#'   entry < f*A reaches the full horizon; the remainder are censored uniformly
#'   on [horizon - (1-f)*A, horizon).
#'
#'   NOTE: an earlier version put the cutoff at A + horizon, which makes
#'   min(horizon, A + horizon - entry) identically equal to horizon for every
#'   entry in [0, A] -- i.e. the mechanism was INERT and censored nobody. The
#'   symptom was Gate 1.6 reporting the same exit-time pile-up (97.2% vs 97.3%)
#'   with accrual on and off.
#' @param scale multiplier on the baseline hazard, set by calibration
#' @param exact_pc optional list(breaks, levels) enabling the closed-form
#'   treated-arm inverse when VE_id is piecewise-constant
simulate_trial <- function(n_c, n_v, h0, ve_id_fn, horizon,
                           theta = 0, ltfu_prob = 0.04, accrual_days = 0,
                           scale = 1, exact_pc = NULL,
                           admin_complete_frac = 0.85) {
  n <- n_c + n_v
  arm <- rep(c(0L, 1L), times = c(n_c, n_v))

  frailty <- if (theta > 0) rgamma(n, shape = 1 / theta, scale = theta) else rep(1, n)

  t_event <- numeric(n)
  t_event[arm == 0L] <- sample_event_times(n_c, h0, NULL, frailty[arm == 0L],
                                           scale = scale)
  t_event[arm == 1L] <- sample_event_times(n_v, h0, ve_id_fn, frailty[arm == 1L],
                                           scale = scale, exact_pc = exact_pc)

  # Loss to follow-up: exponential with rate chosen so P(LTFU <= horizon) = ltfu_prob
  t_ltfu <- if (ltfu_prob > 0) {
    rexp(n, rate = -log1p(-ltfu_prob) / horizon)
  } else {
    rep(Inf, n)
  }

  # Administrative censoring from staggered entry against a fixed data cutoff.
  t_admin <- if (accrual_days > 0) {
    entry <- runif(n, 0, accrual_days)
    pmin(horizon, horizon + admin_complete_frac * accrual_days - entry)
  } else {
    rep(horizon, n)
  }

  t_cens <- pmin(t_ltfu, t_admin)
  time <- pmin(t_event, t_cens)
  status <- as.integer(t_event <= t_cens)

  data.frame(id = seq_len(n), arm = arm, time = time, status = status,
             frailty = frailty, t_event = t_event, t_cens = t_cens)
}

# --- calibration ------------------------------------------------------------
#
# Solve the baseline scale so that the EXPECTED number of observed events over
# both arms matches a target. Expectation is computed by numeric integration
# against the marginal survival and censoring distributions, not by simulation,
# so calibration itself contributes no Monte Carlo noise.

expected_events <- function(scale, n_c, n_v, h0, ve_id_fn, horizon, theta,
                            ltfu_prob, accrual_days, step = 0.25,
                            admin_complete_frac = 0.85) {
  grid <- seq(0, horizon, by = step)
  h0s <- list(cuts = h0$cuts, vals = h0$vals * scale)

  # censoring survival: P(not censored by t)
  ltfu_rate <- if (ltfu_prob > 0) -log1p(-ltfu_prob) / horizon else 0
  s_ltfu <- exp(-ltfu_rate * grid)
  # entry ~ U(0, A), admin follow-up = min(horizon, horizon + f*A - entry), so
  # P(admin > t) = P(entry < horizon + f*A - t) = min(1, (horizon + f*A - t)/A).
  # At t = horizon this equals f, the fraction with complete follow-up.
  s_admin <- if (accrual_days > 0) {
    pmin(1, pmax(0, (horizon + admin_complete_frac * accrual_days - grid) / accrual_days))
  } else {
    rep(1, length(grid))
  }
  s_cens <- s_ltfu * s_admin

  arm_events <- function(n, cumhaz, haz) {
    # marginal (frailty-averaged) survival and hazard
    if (theta > 0) {
      s <- (1 + theta * cumhaz)^(-1 / theta)
      lam <- haz / (1 + theta * cumhaz)
    } else {
      s <- exp(-cumhaz); lam <- haz
    }
    f <- s * lam * s_cens
    n * sum((head(f, -1) + tail(f, -1)) / 2 * step)
  }

  H_C <- pc_integral(h0s, grid); h_C <- pc_eval(h0s, grid)
  H_V <- cumhaz_treated(grid, h0s, ve_id_fn, step = step)
  h_V <- h_C * (1 - ve_id_fn(grid))

  arm_events(n_c, H_C, h_C) + arm_events(n_v, H_V, h_V)
}

calibrate_scale <- function(target_events, n_c, n_v, h0, ve_id_fn, horizon,
                            theta = 0, ltfu_prob = 0.04, accrual_days = 0,
                            interval = c(1e-3, 1e3)) {
  f <- function(s) expected_events(s, n_c, n_v, h0, ve_id_fn, horizon, theta,
                                   ltfu_prob, accrual_days) - target_events
  uniroot(f, interval = interval, tol = 1e-10)$root
}

# --- verified baseline hazard ----------------------------------------------
#
# Inverted from the MATISSE final-analysis placebo cumulative incidences
# (ClinicalTrials.gov NCT04424316) via Lambda(t) = -log(1 - F(t)).
# See notes/source-parameters.md. Resolution below 90 days is not recoverable,
# so 0-90 is a single block -- declared, not smoothed over.

h0_matisse <- list(
  cuts = c(0, 90, 120, 150, 180, 210, 240, 270),
  vals = c(1.905e-4, 2.724e-4, 2.057e-4, 2.069e-4,
           1.042e-4, 0.695e-4, 0.696e-4, 0.817e-4)
)

h0_flat <- list(cuts = c(0), vals = c(2.0e-4))

h0_increasing <- list(
  cuts = c(0, 60, 120, 180, 240, 300),
  vals = c(1.0e-4, 1.5e-4, 2.0e-4, 2.5e-4, 3.0e-4, 3.5e-4)
)
