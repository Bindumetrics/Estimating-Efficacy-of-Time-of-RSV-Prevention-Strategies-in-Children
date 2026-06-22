# tests/test_m1b.R
# Lightweight M1b checks. Run from project root: Rscript tests/test_m1b.R

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))
source(file.path("R", "methods", "m1b.R"))

config <- make_sim_config(m4_interval_width = 15L)
calibrated <- calibrate_all_scenarios(config)
scenario <- calibrated[["A"]]
grid <- seq.int(config$evaluation$grid_start, scenario$spec$horizon)

dat <- simulate_trial(n = 5000L, calibrated = scenario, loss_to_followup = 0, seed = 303)
fit <- fit_m1b(dat, grid = grid, baseline_df = 4L, tvc_df = 3L)
stopifnot(isTRUE(fit$converged))
stopifnot(identical(fit$method, "M1b_supplementary"))
stopifnot(isTRUE(all(fit$estimates$supplementary)))
stopifnot(nrow(fit$estimates) == length(grid))
stopifnot(all(is.finite(fit$estimates$ve)))
stopifnot(all(fit$estimates$lower <= fit$estimates$upper))
stopifnot(all(fit$estimates$day == grid))

few_events <- dat
few_events$status <- 0L
failed <- fit_m1b(few_events, grid = grid, baseline_df = 4L, tvc_df = 3L)
stopifnot(!isTRUE(failed$converged))
stopifnot(grepl("Too few events", failed$failure_reason))

cat("All M1b checks passed.\n")

