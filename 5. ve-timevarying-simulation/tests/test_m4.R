# tests/test_m4.R
# Lightweight M4 checks. Run from project root: Rscript tests/test_m4.R
# Requires a working cmdstanr + CmdStan installation.

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))
source(file.path("R", "methods", "m4.R"))

config <- make_sim_config(m4_interval_width = 15L)
calibrated <- calibrate_all_scenarios(config)
scenario <- calibrated[["B"]]
horizon <- scenario$spec$horizon
grid <- seq.int(config$evaluation$grid_start, horizon)

# --- Binning checks (no Stan needed) ----------------------------------------
dat <- simulate_trial(n = 2000L, calibrated = scenario, loss_to_followup = 0, seed = 808)
ig <- m4_interval_grid(horizon, config$m4$interval_width)
bin <- m4_bin_counts(dat, ig)

stopifnot(nrow(bin) == nrow(ig))
# At-risk counts must be non-increasing across intervals within an arm.
stopifnot(all(diff(bin$n_p) <= 0), all(diff(bin$n_v) <= 0))
# First-interval at-risk equals the arm sizes (no left truncation, start = 0).
stopifnot(bin$n_p[1] == sum(dat$z == 0L), bin$n_v[1] == sum(dat$z == 1L))
# Total binned events equal observed events per arm (no event past horizon).
stopifnot(sum(bin$x_p) == sum(dat$status[dat$z == 0L]))
stopifnot(sum(bin$x_v) == sum(dat$status[dat$z == 1L]))
cat("M4 binning checks passed.\n")

# --- Stan fit checks (skipped if CmdStan unavailable) ------------------------
has_cmdstan <- requireNamespace("cmdstanr", quietly = TRUE) &&
  !is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL))

if (!has_cmdstan) {
  cat("CmdStan not available; skipping M4 sampling checks.\n")
} else {
  model <- m4_stan_model()
  fit <- fit_m4(dat, grid = grid, horizon = horizon,
                interval_width = config$m4$interval_width, model = model,
                seed = 808, iter_warmup = 500L, iter_sampling = 500L,
                chains = 2L, parallel_chains = 2L)
  stopifnot(isTRUE(fit$converged))
  stopifnot(fit$selected_form %in% c("exponential", "erlang3", "power_law"))
  stopifnot(nrow(fit$estimates) == length(grid))
  stopifnot(all(is.finite(fit$estimates$ve)))
  stopifnot(all(fit$estimates$lower <= fit$estimates$upper))
  stopifnot(all(fit$estimates$ve >= -0.5 & fit$estimates$ve <= 1))
  stopifnot(all(is.finite(fit$looic_table$looic)))
  cat("All M4 Stan checks passed. Selected form:", fit$selected_form, "\n")
}
