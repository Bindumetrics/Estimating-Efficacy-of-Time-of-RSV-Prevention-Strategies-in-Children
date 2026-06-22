# tests/test_m3.R
# Lightweight M3 checks. Run from project root: Rscript tests/test_m3.R

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))
source(file.path("R", "methods", "m3.R"))

config <- make_sim_config(m4_interval_width = 15L)
calibrated <- calibrate_all_scenarios(config)
scenario <- calibrated[["A"]]
grid <- seq.int(config$evaluation$grid_start, scenario$spec$horizon)

dat <- simulate_trial(n = 7000L, calibrated = scenario, loss_to_followup = 0, seed = 707)
fit <- fit_m3(dat, grid = grid, bandwidth_days = config$m3$bandwidth_days, kernel = config$m3$kernel)
stopifnot(isTRUE(fit$converged))
stopifnot(identical(fit$method, "M3"))
stopifnot(nrow(fit$estimates) == length(grid))
stopifnot(all(is.finite(fit$estimates$ve)))
stopifnot(all(fit$estimates$lower <= fit$estimates$upper))
stopifnot(all(fit$estimates$day == grid))
stopifnot(all(fit$estimates$local_converged))

few_events <- dat
few_events$status <- 0L
failed <- fit_m3(few_events, grid = grid, bandwidth_days = config$m3$bandwidth_days, kernel = config$m3$kernel)
stopifnot(!isTRUE(failed$converged))
stopifnot(grepl("Fewer than six events", failed$failure_reason))

cat("All M3 checks passed.\n")
