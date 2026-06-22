# methods/m4.R
# -----------------------------------------------------------------------------
# M4: Bayesian parametric model, binomial interval-count likelihood (Section 4.6).
#
# Pipeline:
#   1. Bin the simulated IPD into the fixed interval grid (per-interval at-risk
#      and event counts per arm).
#   2. Fit each of the three candidate VE forms (exponential, Erlang-3,
#      power-law) with Stan via the shared stan/m4_binomial.stan model.
#   3. Compute LOOIC for each form and select the minimum (mirrors the thesis,
#      where the functional form is chosen by LOOIC).
#   4. Return per-form posterior VE(t) curves on the day grid AND the
#      LOOIC-selected curve, with 95% credible intervals.
#
# M4 estimates a risk-scale VE: pi_v,k = pi_p,k * (1 - VE_r,k), where the pi are
# interval-conditional event probabilities. Under rare events this is close to
# (but not identical to) the hazard-scale VE_h(t); the evaluation stage scores M4
# against the DGM-implied interval-conditional VE_r so a pure scale gap is not
# mistaken for method bias.
# -----------------------------------------------------------------------------

# Default weakly informative prior hyperparameters (primary tier).
# These are NOT centred on the simulated truth. The matched Beta(8, 2) prior is
# a sensitivity-tier override, not the default here.
m4_default_priors <- function() {
  list(
    ve0_a = 2, ve0_b = 2,        # VE0 ~ Beta(2, 2)
    w1_sd = 0.05,                # rho / rate ~ Half-Normal(0, 0.05) on day^-1
    w2_meanlog = 0, w2_sdlog = 0.5  # power gamma ~ LogNormal(0, 0.5), median 1
  )
}

m4_form_codes <- function() c(exponential = 1L, erlang3 = 2L, power_law = 3L)

# ---- Interval binning -------------------------------------------------------

# Bin IPD into per-interval, per-arm binomial counts.
# n_*[k] = subjects still event-free and under observation at the start of
# interval k (time >= start_k); x_*[k] = events observed within [start_k, end_k).
m4_bin_counts <- function(data, interval_grid) {
  required <- c("time", "status", "z")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) stop("m4_bin_counts: missing columns: ", paste(missing, collapse = ", "))

  K <- nrow(interval_grid)
  out <- data.frame(
    interval = interval_grid$interval,
    start = interval_grid$start,
    end = interval_grid$end,
    t_mid = (interval_grid$start + interval_grid$end) / 2
  )

  count_arm <- function(arm_value) {
    d <- data[data$z == arm_value, , drop = FALSE]
    n_k <- integer(K)
    x_k <- integer(K)
    for (k in seq_len(K)) {
      a <- interval_grid$start[k]
      b <- interval_grid$end[k]
      at_risk <- d$time >= a
      # Final interval is closed at the top so an event exactly at the horizon
      # (should not occur with administrative censoring) is not dropped.
      in_interval <- if (k == K) (d$time >= a & d$time <= b) else (d$time >= a & d$time < b)
      n_k[k] <- sum(at_risk)
      x_k[k] <- sum(in_interval & d$status == 1L)
    }
    list(n = n_k, x = x_k)
  }

  placebo <- count_arm(0L)
  vaccine <- count_arm(1L)
  out$n_p <- placebo$n
  out$x_p <- placebo$x
  out$n_v <- vaccine$n
  out$x_v <- vaccine$x
  out
}

# ---- Stan model handle (compile once, cache) --------------------------------

# Ensure the C++ toolchain (Rtools on Windows) is visible for model compilation.
# Harmless no-op on platforms / sessions where the toolchain is already on PATH.
m4_ensure_toolchain <- function() {
  if (.Platform$OS.type != "windows") return(invisible())
  candidates <- c("C:/rtools45/usr/bin", "C:/rtools45/x86_64-w64-mingw32.static.posix/bin")
  present <- candidates[dir.exists(candidates)]
  if (length(present) == 0) return(invisible())
  path <- Sys.getenv("PATH")
  add <- present[!vapply(present, function(p) grepl(p, path, fixed = TRUE), logical(1))]
  if (length(add) > 0) Sys.setenv(PATH = paste(paste(add, collapse = ";"), path, sep = ";"))
  invisible()
}

m4_stan_model <- function(stan_file = file.path("stan", "m4_binomial.stan")) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("cmdstanr is required for M4. Install cmdstanr and run cmdstanr::install_cmdstan().")
  }
  m4_ensure_toolchain()
  cmdstanr::cmdstan_model(stan_file)
}

# ---- Fit one candidate form -------------------------------------------------

fit_m4_one_form <- function(model, bin, form, grid, priors,
                            chains = 4L, parallel_chains = 4L,
                            iter_warmup = 1000L, iter_sampling = 1000L,
                            adapt_delta = 0.95, seed = 1L, refresh = 0L) {
  form_code <- m4_form_codes()[[form]]
  stan_data <- list(
    K = nrow(bin),
    n_p = as.integer(bin$n_p), x_p = as.integer(bin$x_p),
    n_v = as.integer(bin$n_v), x_v = as.integer(bin$x_v),
    t_mid = as.numeric(bin$t_mid),
    form = as.integer(form_code),
    G = length(grid),
    t_grid = as.numeric(grid),
    ve0_a = priors$ve0_a, ve0_b = priors$ve0_b,
    w1_sd = priors$w1_sd,
    w2_meanlog = priors$w2_meanlog, w2_sdlog = priors$w2_sdlog
  )

  fit <- tryCatch(
    model$sample(
      data = stan_data,
      chains = chains, parallel_chains = parallel_chains,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      adapt_delta = adapt_delta, seed = seed, refresh = refresh,
      show_messages = FALSE
    ),
    error = function(e) structure(list(error = conditionMessage(e)), class = "m4_form_error")
  )
  if (inherits(fit, "m4_form_error")) {
    return(list(form = form, converged = FALSE, failure_reason = fit$error,
                looic = NA_real_, estimates = NULL, diagnostics = NULL))
  }

  diag <- m4_diagnostics(fit)
  loo_obj <- tryCatch(fit$loo(variables = "log_lik"), error = function(e) NULL)
  looic <- if (is.null(loo_obj)) NA_real_ else loo_obj$estimates["looic", "Estimate"]

  estimates <- m4_summarise_curve(fit, grid, form)

  list(
    form = form,
    converged = diag$ok,
    failure_reason = if (diag$ok) NA_character_ else diag$reason,
    looic = looic,
    loo = loo_obj,
    estimates = estimates,
    diagnostics = diag
  )
}

m4_diagnostics <- function(fit) {
  sm <- tryCatch(fit$diagnostic_summary(quiet = TRUE), error = function(e) NULL)
  draws_summary <- tryCatch(
    fit$summary(variables = c("ve0", "w1"), "rhat", "ess_bulk"),
    error = function(e) NULL
  )
  num_divergent <- if (is.null(sm)) NA_integer_ else sum(sm$num_divergent)
  max_rhat <- if (is.null(draws_summary)) NA_real_ else max(draws_summary$rhat, na.rm = TRUE)
  ok <- isTRUE(is.finite(num_divergent) && num_divergent == 0L &&
                 is.finite(max_rhat) && max_rhat < 1.05)
  reason <- if (ok) NA_character_ else
    sprintf("Convergence flag: divergences=%s, max_rhat=%.3f", num_divergent, max_rhat)
  list(ok = ok, reason = reason, num_divergent = num_divergent, max_rhat = max_rhat)
}

# Posterior summary of VE(t) on the day grid (mean + equal-tailed 95% CrI).
m4_summarise_curve <- function(fit, grid, form, conf_level = 0.95) {
  alpha <- 1 - conf_level
  ve_draws <- tryCatch(
    posterior::as_draws_matrix(fit$draws("ve_grid")),
    error = function(e) NULL
  )
  if (is.null(ve_draws)) return(NULL)

  mean_ve <- apply(ve_draws, 2, mean)
  lower <- apply(ve_draws, 2, stats::quantile, probs = alpha / 2, names = FALSE)
  upper <- apply(ve_draws, 2, stats::quantile, probs = 1 - alpha / 2, names = FALSE)

  data.frame(
    method = "M4",
    form = form,
    day = grid,
    ve = as.numeric(mean_ve),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    conf_level = conf_level
  )
}

# ---- Top-level M4 fit: all three forms + LOOIC selection --------------------

fit_m4 <- function(data, grid, horizon, interval_width = 15L,
                   model = NULL, priors = m4_default_priors(),
                   forms = c("exponential", "erlang3", "power_law"),
                   seed = 1L, ...) {
  if (is.null(model)) model <- m4_stan_model()

  interval_grid <- m4_interval_grid(horizon, interval_width)
  bin <- m4_bin_counts(data, interval_grid)

  if (sum(bin$x_v) + sum(bin$x_p) < 2L) {
    return(m4_failure("Fewer than two events across both arms; M4 is not estimable.", bin))
  }

  per_form <- lapply(forms, function(f) {
    fit_m4_one_form(model, bin, form = f, grid = grid, priors = priors, seed = seed, ...)
  })
  names(per_form) <- forms

  looic_table <- data.frame(
    form = forms,
    looic = vapply(per_form, function(x) x$looic, numeric(1)),
    converged = vapply(per_form, function(x) isTRUE(x$converged), logical(1)),
    stringsAsFactors = FALSE
  )

  eligible <- looic_table[looic_table$converged & is.finite(looic_table$looic), ]
  if (nrow(eligible) == 0L) {
    return(m4_failure("No M4 candidate form converged with a finite LOOIC.", bin,
                      per_form = per_form, looic_table = looic_table))
  }
  selected_form <- eligible$form[which.min(eligible$looic)]
  selected_estimates <- per_form[[selected_form]]$estimates
  if (!is.null(selected_estimates)) selected_estimates$selected <- TRUE

  all_estimates <- do.call(rbind, lapply(per_form, function(x) x$estimates))

  list(
    method = "M4",
    label = "M4, Bayesian binomial interval-count model",
    converged = TRUE,
    failure_reason = NA_character_,
    interval_width = as.integer(interval_width),
    bin = bin,
    per_form = per_form,
    looic_table = looic_table,
    selected_form = selected_form,
    estimates = selected_estimates,    # LOOIC-selected curve (primary)
    all_form_estimates = all_estimates # every form's curve (sensitivity)
  )
}

predict_m4 <- function(fit_object, grid = NULL, conf_level = 0.95) {
  if (!isTRUE(fit_object$converged)) return(fit_object)
  fit_object$estimates
}

m4_failure <- function(reason, bin = NULL, per_form = NULL, looic_table = NULL) {
  list(
    method = "M4",
    label = "M4, Bayesian binomial interval-count model",
    converged = FALSE,
    failure_reason = reason,
    interval_width = NA_integer_,
    bin = bin,
    per_form = per_form,
    looic_table = looic_table,
    selected_form = NA_character_,
    estimates = NULL,
    all_form_estimates = NULL
  )
}
