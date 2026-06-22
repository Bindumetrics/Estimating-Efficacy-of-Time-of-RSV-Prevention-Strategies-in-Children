# dgm.R
# -----------------------------------------------------------------------------
# Data-generating mechanism for the time-varying VE simulation study.
#
# The truth is defined on the hazard scale: h_v(t) = h_0(t) * (1 - VE_h(t)).
# Risk-scale VE_r(t) is derived from the implied cumulative risks and retained
# for fair evaluation of the interval-count Bayesian model (M4).
# -----------------------------------------------------------------------------

source_config_if_needed <- function() {
  if (!exists("make_sim_config", mode = "function")) {
    source(file.path(dirname(sys.frame(1)$ofile %||% getwd()), "config.R"))
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x

clamp <- function(x, lower, upper) pmin(pmax(x, lower), upper)

logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

calibrate_weibull_lambda <- function(target_cumulative_incidence, horizon, shape_k = 1.3) {
  stopifnot(target_cumulative_incidence > 0, target_cumulative_incidence < 1)
  -log(1 - target_cumulative_incidence) / (horizon ^ shape_k)
}

baseline_hazard <- function(t, baseline) {
  lambda <- baseline$lambda
  k <- baseline$shape_k
  ifelse(t <= 0, 0, lambda * k * (t ^ (k - 1)))
}

baseline_cumhaz <- function(t, baseline) {
  lambda <- baseline$lambda
  k <- baseline$shape_k
  ifelse(t <= 0, 0, lambda * (t ^ k))
}

truth_ve_h <- function(t, truth) {
  family <- truth$family
  pars <- truth$parameters
  out <- switch(
    family,
    constant = rep(pars$ve, length(t)),
    exponential = pars$ve0 * exp(-pars$rho * t),
    weibull = pars$ve0 * exp(-((pars$rho * t) ^ pars$gamma)),
    erlang3 = {
      x <- pars$rate * t
      pars$ve0 * exp(-x) * (1 + x + 0.5 * x^2)
    },
    delayed_exponential = {
      elapsed <- pmax(0, t - pars$delay)
      pars$ve0 * exp(-pars$rho * elapsed)
    },
    stop("Unknown VE truth family: ", family)
  )
  clamp(out, 0, 0.999999)
}

vaccine_cumhaz <- function(t, baseline, truth, rel_tol = 1e-8) {
  vapply(t, function(tt) {
    if (tt <= 0) return(0)
    integrate(
      f = function(u) baseline_hazard(u, baseline) * (1 - truth_ve_h(u, truth)),
      lower = 0,
      upper = tt,
      rel.tol = rel_tol,
      subdivisions = 500L
    )$value
  }, numeric(1))
}

cumulative_risks <- function(t, baseline, truth) {
  hc <- baseline_cumhaz(t, baseline)
  hv <- vaccine_cumhaz(t, baseline, truth)
  data.frame(
    day = t,
    placebo_risk = 1 - exp(-hc),
    vaccine_risk = 1 - exp(-hv)
  )
}

truth_ve_r <- function(t, baseline, truth) {
  risks <- cumulative_risks(t, baseline, truth)
  ifelse(risks$placebo_risk <= 0, NA_real_, 1 - risks$vaccine_risk / risks$placebo_risk)
}

anchor_loss <- function(par_unconstrained, spec, baseline) {
  truth <- build_truth_from_unconstrained(par_unconstrained, spec)
  days <- as.numeric(names(spec$anchors))
  target <- as.numeric(spec$anchors)
  implied <- truth_ve_r(days, baseline, truth)
  sum((implied - target)^2) + attr(truth, "penalty")
}

build_truth_from_unconstrained <- function(par, spec) {
  penalty <- 0
  family <- spec$family

  truth <- switch(
    family,
    constant = list(
      family = family,
      parameters = list(ve = spec$fixed$ve)
    ),
    exponential = {
      ve0 <- 0.60 + 0.399 * inv_logit(par[1])
      rho <- exp(par[2])
      list(family = family, parameters = list(ve0 = ve0, rho = rho))
    },
    weibull = {
      ve0 <- 0.60 + 0.399 * inv_logit(par[1])
      rho <- exp(par[2])
      list(family = family, parameters = list(ve0 = ve0, rho = rho, gamma = spec$fixed$gamma))
    },
    erlang3 = {
      ve0 <- 0.60 + 0.399 * inv_logit(par[1])
      rate <- exp(par[2])
      list(family = family, parameters = list(ve0 = ve0, rate = rate))
    },
    delayed_exponential = {
      ve0 <- 0.60 + 0.399 * inv_logit(par[1])
      rho <- exp(par[2])
      list(family = family, parameters = list(ve0 = ve0, rho = rho, delay = spec$fixed$delay))
    },
    stop("Unknown scenario family: ", family)
  )

  attr(truth, "penalty") <- penalty
  truth
}

initial_unconstrained <- function(spec) {
  family <- spec$family
  switch(
    family,
    constant = numeric(0),
    exponential = c(logit((spec$initial$ve0 - 0.60) / 0.399), log(spec$initial$rho)),
    weibull = c(logit((spec$initial$ve0 - 0.60) / 0.399), log(spec$initial$rho)),
    erlang3 = c(logit((spec$initial$ve0 - 0.60) / 0.399), log(spec$initial$rate)),
    delayed_exponential = c(logit((spec$initial$ve0 - 0.60) / 0.399), log(spec$initial$rho)),
    stop("Unknown scenario family: ", family)
  )
}

calibrate_scenario <- function(spec, baseline_shape_k = 1.3) {
  baseline <- list(
    shape_k = baseline_shape_k,
    lambda = calibrate_weibull_lambda(
      target_cumulative_incidence = spec$target_placebo_cumulative_incidence,
      horizon = spec$horizon,
      shape_k = baseline_shape_k
    )
  )

  if (spec$family == "constant") {
    truth <- build_truth_from_unconstrained(numeric(0), spec)
    objective <- 0
    convergence <- 0L
  } else {
    fit <- optim(
      par = initial_unconstrained(spec),
      fn = anchor_loss,
      spec = spec,
      baseline = baseline,
      method = "Nelder-Mead",
      control = list(maxit = 2000, reltol = 1e-12)
    )
    truth <- build_truth_from_unconstrained(fit$par, spec)
    objective <- fit$value
    convergence <- fit$convergence
  }

  list(
    spec = spec,
    baseline = baseline,
    truth = truth,
    calibration = data.frame(
      day = as.numeric(names(spec$anchors)),
      target_cumulative_ve_r = as.numeric(spec$anchors),
      implied_cumulative_ve_r = truth_ve_r(as.numeric(names(spec$anchors)), baseline, truth),
      row.names = NULL
    ),
    optimisation = list(objective = objective, convergence = convergence)
  )
}

calibrate_all_scenarios <- function(config) {
  lapply(config$scenarios, calibrate_scenario, baseline_shape_k = config$baseline$shape_k)
}

make_truth_grid <- function(calibrated, scenario_id, config) {
  horizon <- calibrated$spec$horizon
  days <- seq.int(config$evaluation$grid_start, horizon)
  risks <- cumulative_risks(days, calibrated$baseline, calibrated$truth)
  data.frame(
    scenario = scenario_id,
    day = days,
    ve_h = truth_ve_h(days, calibrated$truth),
    ve_r = 1 - risks$vaccine_risk / risks$placebo_risk,
    placebo_risk = risks$placebo_risk,
    vaccine_risk = risks$vaccine_risk,
    placebo_cumhaz = baseline_cumhaz(days, calibrated$baseline),
    vaccine_cumhaz = vaccine_cumhaz(days, calibrated$baseline, calibrated$truth)
  )
}

allocation_vector <- function(n, ratio) {
  n_v <- round(n * ratio[["vaccine"]] / sum(ratio))
  z <- c(rep(1L, n_v), rep(0L, n - n_v))
  sample(z, length(z), replace = FALSE)
}

sample_control_time <- function(n, baseline) {
  u <- runif(n)
  ((-log(u)) / baseline$lambda) ^ (1 / baseline$shape_k)
}

make_vaccine_inverse <- function(baseline, truth, horizon, grid_size = 5000L) {
  grid <- seq(0, horizon, length.out = grid_size)
  hv <- vaccine_cumhaz(grid, baseline, truth)
  list(
    max_hazard = max(hv),
    draw = function(n) {
      target <- -log(runif(n))
      out <- rep(Inf, n)
      inside <- target <= max(hv)
      out[inside] <- approx(x = hv, y = grid, xout = target[inside], ties = "ordered")$y
      out
    }
  )
}

simulate_trial <- function(n, calibrated, loss_to_followup = 0, seed = NULL, grid_size = 5000L) {
  if (!is.null(seed)) set.seed(seed)
  spec <- calibrated$spec
  z <- allocation_vector(n, spec$randomisation_ratio)

  event_time <- numeric(n)
  control_idx <- which(z == 0L)
  vaccine_idx <- which(z == 1L)

  event_time[control_idx] <- sample_control_time(length(control_idx), calibrated$baseline)
  vaccine_inverse <- make_vaccine_inverse(calibrated$baseline, calibrated$truth, spec$horizon, grid_size)
  event_time[vaccine_idx] <- vaccine_inverse$draw(length(vaccine_idx))

  random_ltfu <- runif(n) < loss_to_followup
  ltfu_time <- rep(Inf, n)
  ltfu_time[random_ltfu] <- runif(sum(random_ltfu), min = 0, max = spec$horizon)

  censor_time <- pmin(spec$horizon, ltfu_time)
  time <- pmin(event_time, censor_time)
  status <- as.integer(event_time <= censor_time)

  data.frame(
    id = seq_len(n),
    z = z,
    time = time,
    status = status,
    event_time = event_time,
    censor_time = censor_time,
    random_ltfu = random_ltfu
  )
}
