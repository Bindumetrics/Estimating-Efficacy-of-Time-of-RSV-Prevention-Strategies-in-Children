# ---------------------------------------------------------------------------
# Scenario definition and the per-replicate runner.
#
# PAIRING is the whole justification for the paired analysis in sec 6.4: every
# estimator must see a byte-identical dataset within a replicate. This is not
# assumed -- the dataset is hashed once and each estimator is handed the same
# object, with the hash re-checked after every call to catch in-place mutation.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(digest) })

ESTIMATORS <- c("cox_tdc", "durham", "tian", "pooled_logistic", "bayes")

make_scenario <- function(name, n_c = 3563, n_v = 3585, horizon = 180,
                          target_events = 200, ve_shape = "gradual",
                          ve_fn = NULL, theta = 0, ltfu_prob = 0.04,
                          accrual_days = 120, h0 = NULL,
                          grid_by = 15, pl_interval = 45, bayes_interval = 30,
                          notes = "") {
  if (is.null(h0)) h0 <- h0_matisse
  if (is.null(ve_fn)) ve_fn <- ve_shapes[[ve_shape]]
  sc <- calibrate_scale(target_events, n_c, n_v, h0, ve_fn, horizon,
                        theta, ltfu_prob, accrual_days)
  grid <- seq(grid_by, horizon, by = grid_by)
  h0s <- list(cuts = h0$cuts, vals = h0$vals * sc)
  structure(list(
    name = name, n_c = n_c, n_v = n_v, horizon = horizon,
    target_events = target_events, ve_shape = ve_shape, ve_fn = ve_fn,
    theta = theta, ltfu_prob = ltfu_prob, accrual_days = accrual_days,
    h0 = h0, h0_scaled = h0s, scale = sc, grid = grid,
    pl_interval = pl_interval, bayes_interval = bayes_interval,
    ve_id_true = ve_fn(grid),
    ve_h_true = ve_h_true(grid, h0s, ve_fn, theta),
    notes = notes), class = "ve_scenario")
}

print.ve_scenario <- function(x, ...) {
  cat(sprintf("<scenario %s> n=%d/%d horizon=%d target=%d shape=%s theta=%.2f ltfu=%.2f accrual=%d scale=%.4f\n",
              x$name, x$n_c, x$n_v, x$horizon, x$target_events, x$ve_shape,
              x$theta, x$ltfu_prob, x$accrual_days, x$scale))
  invisible(x)
}

#' Run every estimator on one replicate of one scenario.
#'
#' Returns a list with `estimates` (long data.frame) and `diagnostics`
#' (one row per method). Failures are recorded with their reason; nothing is
#' silently converted to NA and dropped.
run_replicate <- function(scen, rep_id, estimators = ESTIMATORS,
                          keep_data = FALSE, bayes_seed = NULL) {
  d <- simulate_trial(scen$n_c, scen$n_v, scen$h0, scen$ve_fn, scen$horizon,
                      theta = scen$theta, ltfu_prob = scen$ltfu_prob,
                      accrual_days = scen$accrual_days, scale = scen$scale)
  dat <- d[, c("id", "arm", "time", "status")]
  hash0 <- digest(dat, algo = "xxhash64")

  est <- list(); dg <- list()
  for (m in estimators) {
    r <- switch(m,
      cox_tdc         = try(est_cox_tdc(dat, scen$grid), silent = TRUE),
      durham          = try(est_durham(dat, scen$grid), silent = TRUE),
      tian            = try(est_tian(dat, scen$grid), silent = TRUE),
      pooled_logistic = try(est_pooled_logistic(dat, scen$grid,
                                                interval_width = scen$pl_interval,
                                                horizon = scen$horizon), silent = TRUE),
      bayes           = try(est_bayes(dat, scen$grid, horizon = scen$horizon,
                                      interval_width = scen$bayes_interval,
                                      seed = bayes_seed), silent = TRUE)
    )
    # PAIRING ASSERTION: the dataset must be unchanged after every estimator.
    if (!identical(digest(dat, algo = "xxhash64"), hash0))
      stop(sprintf("estimator '%s' mutated the dataset in scenario '%s' rep %d",
                   m, scen$name, rep_id))

    if (inherits(r, "try-error")) {
      r <- fail_result(scen$grid, NA_real_,
                       paste("uncaught error:", conditionMessage(attr(r, "condition"))))
    }
    est[[m]] <- data.frame(scenario = scen$name, rep = rep_id, method = m,
                           t = r$t, ve_hat = r$ve_hat, lo = r$lo, hi = r$hi)
    d_ <- attr(r, "diagnostics")
    dg[[m]] <- data.frame(
      scenario = scen$name, rep = rep_id, method = m,
      converged = attr(r, "converged"), time_sec = attr(r, "time_sec"),
      note = attr(r, "note"),
      n_na = sum(is.na(r$ve_hat)),
      # test of no waning, harmonised across methods (NA where not defined)
      p_waning = switch(m,
        cox_tdc = d_$wald_p_waning %||% NA_real_,
        durham  = d_$zph_p %||% NA_real_,
        bayes   = d_$p_waning %||% NA_real_,
        NA_real_),
      # surfaced as columns because they are reported per replicate, not just
      # summarised: selection FREQUENCIES and bandwidth stability are findings
      # in their own right (brief sec 9)
      bayes_form = if (m == "bayes") (d_$form_selected %||% NA_integer_) else NA_integer_,
      bayes_elpd_margin_z = if (m == "bayes") (d_$elpd_margin_z %||% NA_real_) else NA_real_,
      bayes_pareto_bad = if (m == "bayes") sum(d_$pareto_k_bad %||% NA_integer_, na.rm = TRUE) else NA_integer_,
      tian_bandwidth = if (m == "tian") (d_$bandwidth %||% NA_real_) else NA_real_,
      tian_h_at_boundary = if (m == "tian") isTRUE(d_$bandwidth_at_boundary) else NA,
      durham_extrapolated = if (m == "durham") sum(d_$extrapolated %||% 0) else NA_integer_,
      extra = I(list(d_)),
      stringsAsFactors = FALSE)
  }
  out <- list(estimates = do.call(rbind, est),
              diagnostics = do.call(rbind, dg),
              n_events = sum(dat$status), data_hash = hash0)
  if (keep_data) out$data <- dat
  out
}
