# ---------------------------------------------------------------------------
# 5. Bayesian parametric waning model (CmdStan)
# ---------------------------------------------------------------------------
# Same interface as the frequentist estimators. Candidate waning forms are
# compared by LOO and the selected form's posterior supplies VE(t) with an
# equal-tailed credible interval.
#
# SELECTION FREQUENCIES, not the modal pick, are what get reported (brief sec 9):
# with ~6-12 interval-level cells these criteria are unstable, so every fit
# records which form won and by what margin.

suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })

# Aggregate individual data to (interval x arm) counts and EXACT person-time.
group_survival <- function(d, horizon, interval_width = 30) {
  brk <- seq(0, horizon, by = interval_width)
  if (tail(brk, 1) < horizon) brk <- c(brk, horizon)
  K <- length(brk) - 1L
  out <- lapply(seq_len(K), function(k) {
    lo <- brk[k]; hi <- brk[k + 1]
    expo <- pmax(0, pmin(d$time, hi) - lo)
    ev <- d$status == 1L & d$time > lo & d$time <= hi
    c(y0 = sum(ev & d$arm == 0L), y1 = sum(ev & d$arm == 1L),
      pt0 = sum(expo[d$arm == 0L]), pt1 = sum(expo[d$arm == 1L]))
  })
  m <- do.call(rbind, out)
  list(K = K, tmid = head(brk, -1) + diff(brk) / 2,
       y0 = as.integer(m[, "y0"]), y1 = as.integer(m[, "y1"]),
       pt0 = m[, "pt0"], pt1 = m[, "pt1"], breaks = brk)
}

.stan_cache <- new.env(parent = emptyenv())

# The compiled model is cached PER PROCESS. A CmdStanModel wraps an external
# pointer, so a cached object exported to a future worker would arrive stale and
# either fail or silently misbehave. Keying on Sys.getpid() means each worker
# builds its own handle; the compiled executable is already on disk, so this
# costs a lookup rather than a recompile.
get_waning_model <- function(path = "stan/waning.stan") {
  key <- paste0("mod_", Sys.getpid())
  if (is.null(.stan_cache[[key]])) {
    .stan_cache[[key]] <- cmdstan_model(path, compile = TRUE)
  }
  .stan_cache[[key]]
}

.eta_draws <- function(dr, form, tgrid) {
  a <- dr$a
  ltm <- log(tgrid / 30)
  if (form == 1) return(outer(a, rep(1, length(tgrid))))
  b <- dr[["b_raw[1]"]]
  if (form == 2) return(outer(a, ltm, function(A, L) A) + outer(b, ltm))
  tau <- dr[["tau_raw[1]"]] * 30
  w <- sapply(tgrid, function(t) 1 - exp(-t / tau))   # draws x grid
  matrix(a, nrow = length(a), ncol = length(tgrid)) +
    (matrix(b, nrow = length(b), ncol = length(tgrid)) -
       matrix(a, nrow = length(a), ncol = length(tgrid))) * w
}

# Defaults raised after the Gate 4 pilot: at 2 chains x 600 sampling draws the
# convergence criterion (rhat < 1.01, min ESS > 400, zero divergences) was met in
# only 70-80% of replicates, with min ESS landing at 381-593. 1500 sampling draws
# roughly doubles ESS. This costs wall-clock time and the projection reflects it;
# the alternative -- relaxing the criterion until the existing runs pass -- would
# be tuning the threshold rather than fixing the problem.
est_bayes <- function(d, grid, horizon = NULL, interval_width = 30,
                      forms = 1:3, chains = 2, iter_warmup = 750,
                      iter_sampling = 1500, level = 0.95, seed = NULL) {
  t0 <- proc.time()[["elapsed"]]
  if (is.null(horizon)) horizon <- max(grid)
  gs <- group_survival(d, horizon, interval_width)
  if (sum(gs$y0) + sum(gs$y1) < 5)
    return(fail_result(grid, proc.time()[["elapsed"]] - t0, "fewer than 5 events"))

  mod <- get_waning_model()
  fits <- list(); elpd <- rep(NA_real_, length(forms)); se_elpd <- elpd
  pk_bad <- rep(NA_integer_, length(forms))
  for (i in seq_along(forms)) {
    f <- forms[i]
    dat <- list(K = gs$K, tmid = gs$tmid, y0 = gs$y0, y1 = gs$y1,
                pt0 = gs$pt0, pt1 = gs$pt1, form = f)
    fit <- try(suppressMessages(mod$sample(
      data = dat, chains = chains, parallel_chains = 1,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE,
      seed = if (is.null(seed)) sample.int(.Machine$integer.max, 1) else seed + f)),
      silent = TRUE)
    if (inherits(fit, "try-error")) next
    ll <- try(fit$draws("log_lik", format = "matrix"), silent = TRUE)
    if (inherits(ll, "try-error")) next
    lo <- try(suppressWarnings(loo::loo(ll)), silent = TRUE)
    if (!inherits(lo, "try-error")) {
      elpd[i] <- lo$estimates["elpd_loo", "Estimate"]
      se_elpd[i] <- lo$estimates["elpd_loo", "SE"]
      # Pareto k > 0.7 means the LOO approximation is unreliable for that point.
      # With only 2K cells this happens often; it is recorded, not suppressed,
      # because it is direct evidence for the brief's sec 9 warning that
      # information criteria are unstable on small grouped data.
      pk_bad[i] <- sum(loo::pareto_k_values(lo) > 0.7)
    }
    fits[[as.character(f)]] <- fit
  }
  el <- proc.time()[["elapsed"]] - t0
  if (length(fits) == 0)
    return(fail_result(grid, el, "all Stan fits failed"))
  if (all(is.na(elpd)))
    return(fail_result(grid, el, "LOO failed for every candidate form"))

  best_i <- which.max(replace(elpd, is.na(elpd), -Inf))
  best <- forms[best_i]
  fit <- fits[[as.character(best)]]

  # convergence diagnostics on the selected fit
  sm <- fit$summary(c("a"))
  rhat_max <- suppressWarnings(max(fit$summary()$rhat, na.rm = TRUE))
  ess_min  <- suppressWarnings(min(fit$summary()$ess_bulk, na.rm = TRUE))
  ndiv <- sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
  converged <- is.finite(rhat_max) && rhat_max < 1.01 && ess_min > 400 && ndiv == 0

  dr <- as_draws_df(fit$draws(c("a",
                                if (best != 1) "b_raw[1]",
                                if (best == 3) "tau_raw[1]")))
  eta <- .eta_draws(dr, best, pmax(grid, 1e-6))
  ve  <- 1 - exp(eta)
  al <- (1 - level) / 2
  qs <- apply(ve, 2, quantile, probs = c(al, 0.5, 1 - al), names = FALSE)

  # Test of no waning -- SELECTION-INDEPENDENT.
  #
  # Earlier this conditioned on the selected form, and set p_waning = 0 whenever
  # the constant form was selected. That was doubly wrong: (i) selecting the
  # constant form is evidence FOR no waning, so it should give a LARGE p-value,
  # not zero, and (ii) making the test depend on a model-selection step that is
  # itself a coin flip (elpd margin ~0.2 SE) inflates the error rate arbitrarily.
  # The observed consequence was an 83% false-rejection rate under constant VE.
  #
  # The fix tests the slope of the LOG-T form (form 2) unconditionally, whichever
  # form is selected for the point estimate. This is the exact Bayesian analogue
  # of the Cox-TDC Wald test on its log-t interaction, so the two are directly
  # comparable, and it does not depend on selection at all. p_waning is the
  # two-sided posterior tail probability that the slope is zero.
  form2_fit <- fits[["2"]]
  p_waning <- if (is.null(form2_fit)) NA_real_ else {
    b <- as_draws_df(form2_fit$draws("b_raw[1]"))[["b_raw[1]"]]
    2 * min(mean(b > 0), mean(b < 0))
  }

  new_result(grid, qs[2, ], qs[1, ], qs[3, ], converged, el,
             if (converged) "" else
               sprintf("rhat=%.3f ess=%.0f divergent=%d", rhat_max, ess_min, ndiv),
             list(form_selected = best, elpd = elpd, se_elpd = se_elpd,
                  elpd_margin = if (sum(!is.na(elpd)) > 1)
                    sort(elpd, decreasing = TRUE)[1] - sort(elpd, decreasing = TRUE)[2]
                    else NA_real_,
                  # margin relative to its own SE: < 1 means the selection is
                  # not distinguishable from a coin flip
                  elpd_margin_z = if (sum(!is.na(elpd)) > 1)
                    (sort(elpd, decreasing = TRUE)[1] - sort(elpd, decreasing = TRUE)[2]) /
                      max(se_elpd, na.rm = TRUE) else NA_real_,
                  pareto_k_bad = pk_bad,
                  rhat_max = rhat_max, ess_min = ess_min, n_divergent = ndiv,
                  p_waning = p_waning, n_intervals = gs$K))
}
