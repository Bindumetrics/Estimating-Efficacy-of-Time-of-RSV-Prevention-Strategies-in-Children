# tests/test_m2.R
# Lightweight M2 checks. Run from project root: Rscript tests/test_m2.R

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))
source(file.path("R", "methods", "m2.R"))

config <- make_sim_config(m4_interval_width = 15L)
calibrated <- calibrate_all_scenarios(config)
scenario <- calibrated[["A"]]
grid <- seq.int(config$evaluation$grid_start, scenario$spec$horizon)

dat <- simulate_trial(n = 7000L, calibrated = scenario, loss_to_followup = 0, seed = 606)
fit <- fit_m2(dat, grid = grid, basis_dimension = config$m2$basis_dimension)
stopifnot(isTRUE(fit$converged))
stopifnot(identical(fit$method, "M2"))
stopifnot(nrow(fit$estimates) == length(grid))
stopifnot(all(is.finite(fit$estimates$ve)))
stopifnot(all(fit$estimates$lower <= fit$estimates$upper))
stopifnot(all(fit$estimates$day == grid))
stopifnot(nrow(fit$event_curve) == sum(dat$status == 1L))

few_events <- dat
few_events$status <- 0L
failed <- fit_m2(few_events, grid = grid, basis_dimension = config$m2$basis_dimension)
stopifnot(!isTRUE(failed$converged))
stopifnot(grepl("Fewer than six events", failed$failure_reason))

cat("All M2 checks passed.\n")
