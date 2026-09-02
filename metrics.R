# ---------------------------------------------------------------------------
# Metrics and analysis (brief sec 6).
#
# Everything is POINTWISE in t first; aggregates are derived from the pointwise
# quantities, never the other way round. Every reported number carries a Monte
# Carlo standard error, and comparisons smaller than their MCSE are labelled
# inconclusive rather than reported as differences.
#
# MCSE formulas follow Morris, White & Crowther (2019, Stat Med 38:2074-2102):
#   bias      MCSE = sd(theta_hat) / sqrt(R)
#   empSE     MCSE = empSE / sqrt(2(R-1))
#   MSE       MCSE = sd((theta_hat - theta)^2) / sqrt(R)
#   coverage  MCSE = sqrt(p(1-p)/R)
# RMSE and its MCSE are obtained from MSE by the delta method.
# ---------------------------------------------------------------------------

# `estimates` is the long frame from run_replicate(); `truth` is a data.frame
# with columns t and ve_true carrying the ESTIMAND (VE_h, not VE_id).

pointwise_metrics <- function(estimates, truth, blowup_rule = 1.0) {
  x <- merge(estimates, truth, by = "t")
  x$excluded <- is.finite(x$ve_hat) & abs(x$ve_hat) > blowup_rule
  do.call(rbind, lapply(split(x, list(x$scenario, x$method, x$t), drop = TRUE), function(g) {
    keep <- !is.na(g$ve_hat) & !g$excluded
    R <- sum(keep)
    # Always return the SAME columns, whether or not the cell is estimable, so
    # do.call(rbind, ...) can stack every row. A cell with fewer than 2 usable
    # replicates (all methods failed, or a method not run in this scenario) gets
    # NA metrics rather than a short row with a different schema -- the earlier
    # short-circuit returned 5 columns against the full 20 and broke the rbind.
    if (R < 2) {
      e <- cov <- wid <- numeric(0); mse <- p <- emp <- NA_real_
    } else {
      e   <- g$ve_hat[keep] - g$ve_true[keep]
      cov <- as.numeric(g$lo[keep] <= g$ve_true[keep] & g$ve_true[keep] <= g$hi[keep])
      wid <- g$hi[keep] - g$lo[keep]
      mse <- mean(e^2); p <- mean(cov); emp <- sd(g$ve_hat[keep])
    }
    data.frame(
      scenario = g$scenario[1], method = g$method[1], t = g$t[1],
      R_used = R, R_total = nrow(g),
      ve_true = g$ve_true[1],
      bias = if (R < 2) NA_real_ else mean(e),
      bias_mcse = if (R < 2) NA_real_ else sd(e) / sqrt(R),
      bias_median = if (R < 2) NA_real_ else median(e),
      rmse = if (R < 2) NA_real_ else sqrt(mse),
      rmse_mcse = if (R < 2) NA_real_ else sd(e^2) / sqrt(R) / (2 * sqrt(mse)),
      emp_se = emp,
      emp_se_mcse = if (R < 2) NA_real_ else emp / sqrt(2 * (R - 1)),
      coverage = p,
      coverage_mcse = if (R < 2) NA_real_ else sqrt(p * (1 - p) / R),
      width = if (R < 2) NA_real_ else mean(wid),
      width_mcse = if (R < 2) NA_real_ else sd(wid) / sqrt(R),
      width_median = if (R < 2) NA_real_ else median(wid),
      na_rate = mean(is.na(g$ve_hat)),
      blowup_rate = mean(g$excluded, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }))
}

# Aggregate summaries. Both MEAN and MEDIAN are reported: the mean is the
# conventional definition, the median is robust to heavy tails, and a material
# gap between them is itself a finding about replicate-level stability.
summary_metrics <- function(pw) {
  do.call(rbind, lapply(split(pw, list(pw$scenario, pw$method), drop = TRUE), function(g) {
    data.frame(
      scenario = g$scenario[1], method = g$method[1],
      mean_abs_bias = mean(abs(g$bias)), median_abs_bias = median(abs(g$bias)),
      integrated_sq_err = mean(g$rmse^2),
      mean_rmse = mean(g$rmse), median_rmse = median(g$rmse),
      overall_coverage = mean(g$coverage),
      worst_coverage = min(g$coverage),
      worst_coverage_t = g$t[which.min(g$coverage)],
      mean_width = mean(g$width), median_width = median(g$width_median),
      mean_na_rate = mean(g$na_rate), mean_blowup = mean(g$blowup_rate),
      stringsAsFactors = FALSE)
  }))
}

# Paired comparison: per-replicate differences, which is the entire payoff of
# forcing every estimator to see byte-identical data. The MCSE of the paired
# difference is much smaller than either marginal MCSE.
paired_comparisons <- function(estimates, truth, metric = c("abs_err", "sq_err")) {
  metric <- match.arg(metric)
  x <- merge(estimates, truth, by = "t")
  x$loss <- if (metric == "abs_err") abs(x$ve_hat - x$ve_true) else (x$ve_hat - x$ve_true)^2
  # one loss value per (scenario, method, replicate), averaged over the grid
  agg <- aggregate(loss ~ scenario + method + rep, data = x, FUN = mean, na.rm = TRUE)
  out <- list()
  for (sc in unique(agg$scenario)) {
    a <- agg[agg$scenario == sc, ]
    w <- reshape(a, idvar = "rep", timevar = "method", direction = "wide",
                 drop = "scenario")
    ms <- sub("^loss\\.", "", setdiff(names(w), "rep"))
    for (i in seq_along(ms)) for (j in seq_along(ms)) {
      if (i >= j) next
      d <- w[[paste0("loss.", ms[i])]] - w[[paste0("loss.", ms[j])]]
      d <- d[is.finite(d)]
      if (length(d) < 2) next
      mc <- sd(d) / sqrt(length(d))
      out[[length(out) + 1]] <- data.frame(
        scenario = sc, method_a = ms[i], method_b = ms[j], metric = metric,
        n_paired = length(d), mean_diff = mean(d), mcse = mc,
        lo = mean(d) - 1.96 * mc, hi = mean(d) + 1.96 * mc,
        conclusive = abs(mean(d)) > 1.96 * mc,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

# Type I error (constant-VE scenarios) and power (declining scenarios), from
# each method's test of no waning. Reproducing Haber's comparison requires this;
# accuracy metrics alone do not engage with his result.
power_typeI <- function(diagnostics, alpha = 0.05) {
  d <- diagnostics[is.finite(diagnostics$p_waning), ]
  if (!nrow(d)) return(NULL)
  do.call(rbind, lapply(split(d, list(d$scenario, d$method), drop = TRUE), function(g) {
    R <- nrow(g); p <- mean(g$p_waning < alpha)
    data.frame(scenario = g$scenario[1], method = g$method[1], alpha = alpha,
               R = R, reject_rate = p, mcse = sqrt(p * (1 - p) / R),
               stringsAsFactors = FALSE)
  }))
}

# Smoothness / stability of the estimated curves: mean absolute second
# difference across the grid, per replicate, then averaged.
curve_stability <- function(estimates) {
  do.call(rbind, lapply(split(estimates, list(estimates$scenario, estimates$method,
                                              estimates$rep), drop = TRUE), function(g) {
    g <- g[order(g$t), ]
    v <- g$ve_hat
    if (sum(!is.na(v)) < 3) return(NULL)
    d2 <- diff(diff(v))
    data.frame(scenario = g$scenario[1], method = g$method[1], rep = g$rep[1],
               mean_abs_d2 = mean(abs(d2), na.rm = TRUE), stringsAsFactors = FALSE)
  }))
}

# Ranking. The full per-criterion rank matrix is printed so the summary is
# auditable, and the mean-rank column is verified to reproduce from it.
rank_matrix <- function(sm, criteria = c("mean_abs_bias", "mean_rmse",
                                         "mean_width", "cov_dev"),
                        weights = NULL) {
  sm$cov_dev <- abs(sm$overall_coverage - 0.95)
  out <- do.call(rbind, lapply(split(sm, sm$scenario), function(g) {
    r <- sapply(criteria, function(cr) rank(g[[cr]], ties.method = "average"))
    r <- matrix(r, nrow = nrow(g), dimnames = list(g$method, criteria))
    w <- if (is.null(weights)) rep(1 / length(criteria), length(criteria)) else weights
    data.frame(scenario = g$scenario[1], method = g$method, r,
               mean_rank = as.vector(r %*% w), stringsAsFactors = FALSE)
  }))
  # verification: mean_rank must reproduce exactly from the printed columns
  w <- if (is.null(weights)) rep(1 / length(criteria), length(criteria)) else weights
  chk <- as.matrix(out[, criteria]) %*% w
  attr(out, "mean_rank_reproduces") <- isTRUE(all.equal(as.vector(chk), out$mean_rank))
  out
}

# Sensitivity of conclusions to the blow-up exclusion rule. If the ranking flips
# when the threshold moves, the ranking is an artefact of the rule and must be
# reported as such.
blowup_sensitivity <- function(estimates, truth, rules = c(0.5, 1.0, 2.0, Inf)) {
  do.call(rbind, lapply(rules, function(r) {
    pw <- pointwise_metrics(estimates, truth, blowup_rule = r)
    sm <- summary_metrics(pw)
    sm$rule <- r
    sm[, c("rule", "scenario", "method", "mean_abs_bias", "mean_rmse",
           "overall_coverage", "mean_blowup")]
  }))
}
