# tests/test_m1.R
# Lightweight M1 checks. Run from project root: Rscript tests/test_m1.R

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))
source(file.path("R", "methods", "m1.R"))

config <- make_sim_config(m4_interval_width = 15L)
calibrated <- calibrate_all_scenarios(config)
scenario <- calibrated[["A"]]
grid <- seq.int(config$evaluation$grid_start, scenario$spec$horizon)

dat <- simulate_trial(n = 5000L, calibrated = scenario, loss_to_followup = 0, seed = 202)
fit <- fit_m1(dat, grid = grid)
stopifnot(isTRUE(fit$converged))
stopifnot(all(c("z", "log_time_interaction") %in% names(fit$coefficients)))
stopifnot(nrow(fit$estimates) == length(grid))
stopifnot(all(is.finite(fit$estimates$ve)))
stopifnot(all(fit$estimates$lower <= fit$estimates$upper))
stopifnot(all(fit$estimates$day == grid))

few_events <- dat
few_events$status <- 0L
failed <- fit_m1(few_events, grid = grid)
stopifnot(!isTRUE(failed$converged))
stopifnot(grepl("Fewer than two events", failed$failure_reason))

cat("All M1 checks passed.\n")

