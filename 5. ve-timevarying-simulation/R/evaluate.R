# evaluate.R
# -----------------------------------------------------------------------------
# Stage 3: grid harmonisation, scale handling, and estimation metrics.
#
# Every method is scored on the common integer-day grid (1..horizon). The
# decisive fairness point is the SCALE of the truth each method targets:
#
#   * M1, M1b, M2, M3 estimate a hazard-scale VE_h(t)  -> scored against VE_h(t).
#   * M4 estimates a risk-scale VE that enters the binomial likelihood as
#     pi_v,k = pi_p,k (1 - VE_r,k). The matching DGM truth is the
#     interval-conditional risk-scale VE (1 - conditional vaccine risk /
#     conditional placebo risk per interval), NOT the cumulative VE_r and NOT
#     VE_h. Scoring M4 against VE_h would let a pure scale gap masquerade as bias.
#
# Metrics per (method, scenario, n, censoring, grid-day): signed bias, MAE, MSE,
# RMSE, coverage of the 95% interval, mean interval width; ISE over the horizon
# per replication; and convergence/failure rate as its own outcome. For M3 the
# intervals are the pointwise sandwich intervals (a simultaneous band would
# over-cover at any single day).
#
# Two further comparison layers (RQ4 directly):
#   * COMMON-SCALE comparison. To rank all four methods head-to-head we also
#     score every method against ONE common reference curve (default the
#     hazard-scale VE_h(t), the VE the DGM was built from). M4's risk-scale
#     estimate then carries a small, known scale offset; `scale_gap_table()`
#     quantifies max|VE_h - VE_r| per scenario so that offset is transparent and
#     not mistaken for estimation error. The own-scale metrics (each method vs
#     the estimand it targets) remain the primary, bias-pure summary.
#   * M4 PER-FORM scoring. M4 is expanded into the LOOIC-selected curve ("M4")
#     plus one pseudo-method per candidate form ("M4_exponential", "M4_erlang3",
#     "M4_power_law"), so the misspecification tier reports every form against
#     every truth, not only the selected pairing.
# -----------------------------------------------------------------------------

# Hazard- vs risk-scale target for each method.
method_scale <- function(method) {
  ifelse(grepl("^M4", toupper(method)), "risk", "hazard")
}

# Map each integer day to the interval it falls in, using the (start, end]
# convention so day 15 belongs to the first 15-day interval.
day_to_interval <- function(days, interval_grid) {
  idx <- findInterval(days, interval_grid$end, left.open = TRUE) + 1L
  pmin(idx, nrow(interval_grid))
}

# DGM-implied interval-conditional risk-scale VE that M4 targets, per interval,
# mapped onto the day grid as a step function.
truth_ve_r_interval <- function(days, calibrated, interval_grid) {
  d_Hc <- baseline_cumhaz(interval_grid$end, calibrated$baseline) -
    baseline_cumhaz(interval_grid$start, calibrated$baseline)
  d_Hv <- vaccine_cumhaz(interval_grid$end, calibrated$baseline, calibrated$truth) -
    vaccine_cumhaz(interval_grid$start, calibrated$baseline, calibrated$truth)
  pi_p <- 1 - exp(-d_Hc)
  pi_v <- 1 - exp(-d_Hv)
  ve_k <- 1 - pi_v / pi_p
  ve_k[day_to_interval(days, interval_grid)]
}

# Both truth targets on the day grid for one scenario.
truth_targets <- function(calibrated, scenario_id, config) {
  grid <- seq.int(config$evaluation$grid_start, calibrated$spec$horizon)
  tg <- make_truth_grid(calibrated, scenario_id, config)  # supplies ve_h
  interval_grid <- m4_interval_grid(calibrated$spec$horizon, config$m4$interval_width)
  data.frame(
    scenario = scenario_id,
    day = grid,
    truth_hazard = tg$ve_h,
    truth_risk = truth_ve_r_interval(grid, calibrated, interval_grid)
  )
}

# Scale gap between the two truth scales for one scenario: how far the
# hazard-scale VE_h and the interval-conditional risk-scale VE_r diverge. Small
# values justify a common-scale four-way comparison.
scale_gap_summary <- function(calibrated, scenario_id, config) {
  tt <- truth_targets(calibrated, scenario_id, config)
  gap <- abs(tt$truth_hazard - tt$truth_risk)
  data.frame(
    scenario = scenario_id,
    horizon = calibrated$spec$horizon,
    max_abs_gap = max(gap),
    mean_abs_gap = mean(gap),
    gap_at_horizon = gap[length(gap)]
  )
}

# Scale gap for every scenario (calibrates internally).
scale_gap_table <- function(config) {
  calibrated <- calibrate_all_scenarios(config)
  do.call(rbind, lapply(names(calibrated), function(id) {
    scale_gap_summary(calibrated[[id]], id, config)
  }))
}

# ---- Harmonisation ----------------------------------------------------------

# Pull one method's estimates onto the full grid, filling failed fits / days
# with NA. Returns method, day, ve, lower, upper (always length(grid) rows).
extract_estimates <- function(fit, grid) {
  method <- fit$method %||% "unknown"
  base <- data.frame(method = method, day = grid,
                     ve = NA_real_, lower = NA_real_, upper = NA_real_,
                     converged = isTRUE(fit$converged))
  est <- fit$estimates
  if (is.null(est) || nrow(est) == 0) return(base)
  keep <- est[, intersect(c("day", "ve", "lower", "upper"), names(est)), drop = FALSE]
  m <- match(grid, keep$day)
  base$ve <- keep$ve[m]
  base$lower <- keep$lower[m]
  base$upper <- keep$upper[m]
  base
}

harmonise_estimates <- function(fits, grid) {
  out <- do.call(rbind, lapply(fits, extract_estimates, grid = grid))
  out$scale <- method_scale(out$method)
  rownames(out) <- NULL
  out
}

# Expand an M4 fit into pseudo-method fits: the LOOIC-selected curve ("M4") plus
# one per candidate form ("M4_<form>"). Each is a minimal fit-like list that
# harmonise_estimates / extract_estimates can consume. Non-M4 fits pass through.
expand_m4_forms <- function(m4_fit) {
  out <- list(M4 = list(method = "M4", converged = isTRUE(m4_fit$converged),
                        estimates = m4_fit$estimates))
  if (!is.null(m4_fit$per_form)) {
    for (f in names(m4_fit$per_form)) {
      pf <- m4_fit$per_form[[f]]
      label <- paste0("M4_", f)
      est <- pf$estimates
      if (!is.null(est)) est$method <- label
      out[[label]] <- list(method = label, converged = isTRUE(pf$converged), estimates = est)
    }
  }
  out
}

# Replace any M4 fit (that carries per-form results) with its expansion.
expand_fits_for_scoring <- function(fits, score_m4_per_form) {
  if (!isTRUE(score_m4_per_form)) return(fits)
  out <- list()
  for (nm in names(fits)) {
    f <- fits[[nm]]
    if (identical(f$method, "M4") && !is.null(f$per_form)) {
      out <- c(out, expand_m4_forms(f))
    } else {
      out <- c(out, stats::setNames(list(f), nm))
    }
  }
  out
}

# ---- Scoring one replication ------------------------------------------------

# Attach the own-scale truth (each method vs its own estimand) and the common
# reference truth (all methods vs one curve, for the four-way comparison).
attach_truth <- function(harmonised, targets, common_scale = "hazard") {
  tr <- targets[match(harmonised$day, targets$day), ]
  harmonised$truth <- ifelse(harmonised$scale == "risk", tr$truth_risk, tr$truth_hazard)
  harmonised$truth_common <- if (common_scale == "risk") tr$truth_risk else tr$truth_hazard
  harmonised
}

score_estimates <- function(scored) {
  scored$available <- is.finite(scored$ve) & is.finite(scored$truth)
  scored$width <- scored$upper - scored$lower

  # Own-scale metrics: each method against the estimand it targets (bias-pure).
  scored$error <- scored$ve - scored$truth                       # signed bias
  scored$abs_error <- abs(scored$error)
  scored$sq_error <- scored$error^2
  scored$covered <- as.integer(is.finite(scored$lower) & is.finite(scored$upper) &
                                 scored$truth >= scored$lower & scored$truth <= scored$upper)
  scored$covered[!scored$available] <- NA_integer_

  # Common-scale metrics: every method against one reference curve, so all four
  # methods are directly rankable. Carries each method's scale offset.
  scored$error_c <- scored$ve - scored$truth_common
  scored$abs_error_c <- abs(scored$error_c)
  scored$sq_error_c <- scored$error_c^2
  scored$covered_c <- as.integer(is.finite(scored$lower) & is.finite(scored$upper) &
                                   scored$truth_common >= scored$lower & scored$truth_common <= scored$upper)
  scored$covered_c[!scored$available] <- NA_integer_
  scored
}

# Trapezoidal integral with unit (or arbitrary) day spacing.
trapz <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  o <- order(x); x <- x[o]; y <- y[o]
  sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
}

# Integrated squared error per method for one replication, over available days.
replication_ise <- function(scored) {
  methods <- unique(scored$method)
  do.call(rbind, lapply(methods, function(m) {
    s <- scored[scored$method == m & scored$available, ]
    horizon_days <- sum(scored$method == m)
    data.frame(
      method = m,
      scale = method_scale(m),
      ise = trapz(s$day, s$sq_error),
      iae = trapz(s$day, s$abs_error),
      frac_days_available = if (horizon_days > 0) nrow(s) / horizon_days else 0
    )
  }))
}

# Score all methods for one simulated dataset (one replication of one cell).
# common_scale: the single reference curve all methods are also scored against
# ("hazard" = VE_h, the DGM's defining truth; "risk" = interval VE_r).
# score_m4_per_form: also score each M4 candidate form as its own pseudo-method.
evaluate_methods <- function(fits, calibrated, scenario_id, config,
                             replication = NA_integer_, n = NA_integer_,
                             loss_to_followup = NA_real_,
                             common_scale = "hazard", score_m4_per_form = TRUE) {
  grid <- seq.int(config$evaluation$grid_start, calibrated$spec$horizon)
  targets <- truth_targets(calibrated, scenario_id, config)

  fits <- expand_fits_for_scoring(fits, score_m4_per_form)
  scored <- score_estimates(attach_truth(harmonise_estimates(fits, grid), targets,
                                         common_scale = common_scale))
  cell <- data.frame(scenario = scenario_id, n = n, loss_to_followup = loss_to_followup,
                     replication = replication)
  scored <- cbind(cell[rep(1, nrow(scored)), ], scored, row.names = NULL)

  ise <- replication_ise(scored)
  ise <- cbind(cell[rep(1, nrow(ise)), ], ise, row.names = NULL)

  list(scored = scored, ise = ise)
}

# ---- Aggregation across replications ---------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# Per (method, cell, day) metrics, excluding non-available (failed) estimates
# from point summaries but recording how many replications contributed.
aggregate_per_day <- function(scored_stack) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required for aggregation.")
  keys <- c("scenario", "n", "loss_to_followup", "method", "scale", "day")
  dplyr::summarise(
    dplyr::group_by(dplyr::filter(scored_stack, .data$available),
                    dplyr::across(dplyr::all_of(keys))),
    truth = mean(.data$truth),
    mean_ve = mean(.data$ve, na.rm = TRUE),               # mean estimated curve (for plots)
    sd_ve = stats::sd(.data$ve, na.rm = TRUE),
    mean_lower = mean(.data$lower, na.rm = TRUE),         # mean interval limits (for bands)
    mean_upper = mean(.data$upper, na.rm = TRUE),
    signed_bias = mean(.data$error),
    mae = mean(.data$abs_error),
    mse = mean(.data$sq_error),
    rmse = sqrt(mean(.data$sq_error)),
    # Robust companions: mean bias / MSE are dominated by the occasional
    # degenerate small-n fit (e.g. M1's log-time extrapolation blowing up).
    # Median bias and median absolute error are robust to those; n_extreme
    # counts replications whose VE estimate left the plausible range.
    median_bias = stats::median(.data$error),
    median_abs_error = stats::median(.data$abs_error),
    n_extreme = sum(abs(.data$ve) > 1.5, na.rm = TRUE),
    coverage = mean(.data$covered, na.rm = TRUE),
    mean_width = mean(.data$width, na.rm = TRUE),
    # Common-scale metrics (every method vs the same reference curve): the basis
    # for the head-to-head four-way ranking. Includes each method's scale offset.
    truth_common = mean(.data$truth_common),
    signed_bias_common = mean(.data$error_c),
    median_bias_common = stats::median(.data$error_c),
    mse_common = mean(.data$sq_error_c),
    rmse_common = sqrt(mean(.data$sq_error_c)),
    coverage_common = mean(.data$covered_c, na.rm = TRUE),
    n_replications_used = dplyr::n(),
    .groups = "drop"
  )
}

# Convergence / failure rate as its own outcome: per (method, cell), the
# fraction of replications that produced a usable curve (>= 1 finite day), plus
# the average fraction of grid days available (captures partial M3 failures).
aggregate_convergence <- function(scored_stack) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required for aggregation.")
  keys <- c("scenario", "n", "loss_to_followup", "method", "replication")
  per_rep <- dplyr::summarise(
    dplyr::group_by(scored_stack, dplyr::across(dplyr::all_of(keys))),
    any_available = as.integer(any(.data$available)),
    frac_days_available = mean(.data$available),
    .groups = "drop"
  )
  cell_keys <- c("scenario", "n", "loss_to_followup", "method")
  dplyr::summarise(
    dplyr::group_by(per_rep, dplyr::across(dplyr::all_of(cell_keys))),
    n_replications = dplyr::n(),
    convergence_rate = mean(.data$any_available),
    mean_frac_days_available = mean(.data$frac_days_available),
    .groups = "drop"
  )
}

# Mean ISE / IAE per (method, cell) over replications where it is finite.
aggregate_ise <- function(ise_stack) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required for aggregation.")
  keys <- c("scenario", "n", "loss_to_followup", "method", "scale")
  dplyr::summarise(
    dplyr::group_by(ise_stack, dplyr::across(dplyr::all_of(keys))),
    mean_ise = mean(.data$ise, na.rm = TRUE),
    mean_iae = mean(.data$iae, na.rm = TRUE),
    n_ise = sum(is.finite(.data$ise)),
    .groups = "drop"
  )
}

# Pairwise relative efficiency = ratio of horizon-integrated MSE between methods.
#   * own-scale (scale = "hazard"/"risk", metric_col = "mse"): methods that share
#     a truth scale are directly comparable on their bias-pure MSE.
#   * common-scale (scale = NULL, metric_col = "mse_common"): all four methods
#     ranked head-to-head against one reference curve. This is the RQ4 ranking;
#     read it together with scale_gap_table() so M4's scale offset is visible.
# By default the per-form M4 pseudo-methods are dropped from the matrix; set
# include_m4_forms = TRUE to keep them.
relative_efficiency <- function(per_day_agg, scenario_id, n, loss_to_followup,
                                scale = "hazard", metric_col = "mse",
                                include_m4_forms = FALSE) {
  sub <- per_day_agg[per_day_agg$scenario == scenario_id & per_day_agg$n == n &
                       per_day_agg$loss_to_followup == loss_to_followup, ]
  if (!is.null(scale)) sub <- sub[sub$scale == scale, ]
  if (!include_m4_forms) sub <- sub[!grepl("^M4_", sub$method), ]
  if (nrow(sub) == 0) return(NULL)
  methods <- sort(unique(sub$method))
  imse <- vapply(methods, function(m) {
    s <- sub[sub$method == m, ]
    trapz(s$day, s[[metric_col]]) / (max(s$day) - min(s$day))   # mean metric over the horizon
  }, numeric(1))
  ratio <- outer(imse, imse, "/")                      # ratio[i, j] = MSE_i / MSE_j
  dimnames(ratio) <- list(method = methods, reference = methods)
  list(integrated_mse = imse, ratio = ratio, scale = if (is.null(scale)) "common" else scale,
       metric = metric_col)
}

# Convenience: the four-way head-to-head ranking on the common reference scale.
relative_efficiency_common <- function(per_day_agg, scenario_id, n, loss_to_followup,
                                       include_m4_forms = FALSE) {
  relative_efficiency(per_day_agg, scenario_id, n, loss_to_followup,
                      scale = NULL, metric_col = "mse_common",
                      include_m4_forms = include_m4_forms)
}
