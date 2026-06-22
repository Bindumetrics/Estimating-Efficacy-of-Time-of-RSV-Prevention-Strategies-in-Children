# tests/test_stage1_dgm.R
# Lightweight base-R checks for Stage 1. Run from project root:
# Rscript tests/test_stage1_dgm.R

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))

config <- make_sim_config(m4_interval_width = 15L)
calibrated <- calibrate_all_scenarios(config)

for (id in names(calibrated)) {
  x <- calibrated[[id]]
  stopifnot(x$optimisation$convergence == 0L)

  horizon <- x$spec$horizon
  placebo_risk <- 1 - exp(-baseline_cumhaz(horizon, x$baseline))
  stopifnot(abs(placebo_risk - x$spec$target_placebo_cumulative_incidence) < 1e-10)

  days <- seq.int(1L, horizon)
  ve_h <- truth_ve_h(days, x$truth)
  ve_r <- truth_ve_r(days, x$baseline, x$truth)
  stopifnot(all(is.finite(ve_h)), all(ve_h >= 0), all(ve_h < 1))
  stopifnot(all(is.finite(ve_r)), all(ve_r >= 0), all(ve_r < 1))

  anchors <- x$calibration
  tolerance <- if (id == "E") 0.03 else 0.01
  stopifnot(max(abs(anchors$implied_cumulative_ve_r - anchors$target_cumulative_ve_r)) < tolerance)

  dat <- simulate_trial(n = 500, calibrated = x, loss_to_followup = 0.05, seed = 100 + match(id, names(calibrated)))
  stopifnot(nrow(dat) == 500L)
  stopifnot(all(dat$time >= 0), all(dat$time <= horizon))
  stopifnot(all(dat$status %in% c(0L, 1L)))
  stopifnot(all(dat$z %in% c(0L, 1L)))
}

intervals <- m4_interval_grid(180, config$m4$interval_width)
stopifnot(intervals$start[1] == 0L)
stopifnot(tail(intervals$end, 1L) == 180L)
stopifnot(all(intervals$end > intervals$start))

cat("All Stage 1 DGM checks passed.\n")
