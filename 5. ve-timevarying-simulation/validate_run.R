# validate_run.R
# Stage 4 validation: run the full orchestration on a tiny SMOKE design that
# includes all five methods (M4 with reduced Stan sampling). This confirms the
# sweep wiring end-to-end before any full-scale run. It is NOT a results run.
#
# Run from project root:  Rscript validate_run.R

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))
source(file.path("R", "evaluate.R"))
source(file.path("R", "methods", "m1.R"))
source(file.path("R", "methods", "m1b.R"))
source(file.path("R", "methods", "m2.R"))
source(file.path("R", "methods", "m3.R"))
source(file.path("R", "methods", "m4.R"))
source(file.path("R", "run_simulation.R"))

config <- make_sim_config(m4_interval_width = 15L)

# Smoke design: well-specified-for-M4 exponential (B) and Erlang-3 (D), two
# small sample sizes, no loss, a few replications each.
design <- smoke_design(scenarios = c("B", "D"), sample_sizes = c(500L, 2000L),
                       losses = 0, n_replications = 3L)
controls <- smoke_controls()   # all 5 methods, M4 = 2 chains x 300+300

cat("Smoke design:\n"); print(design)

t0 <- Sys.time()
res <- run_simulation(design, config, controls = controls, parallel = FALSE,
                      out_dir = file.path("outputs", "sweep_smoke"),
                      save_results = TRUE, save_raw_scored = TRUE, verbose = TRUE)
cat(sprintf("\nElapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

summary_days <- config$evaluation$summary_days
pd <- res$per_day
pd_sd <- pd[pd$day %in% summary_days, ]

write.csv(pd, file.path("outputs", "sweep_smoke", "smoke_per_day.csv"), row.names = FALSE)
write.csv(res$convergence, file.path("outputs", "sweep_smoke", "smoke_convergence.csv"), row.names = FALSE)
write.csv(res$ise, file.path("outputs", "sweep_smoke", "smoke_ise.csv"), row.names = FALSE)
if (!is.null(res$m4_selection)) {
  write.csv(res$m4_selection, file.path("outputs", "sweep_smoke", "smoke_m4_selection.csv"), row.names = FALSE)
}

cat("\nConvergence rates by cell and method:\n")
print(res$convergence[, c("scenario", "n", "method", "convergence_rate", "mean_frac_days_available")],
      row.names = FALSE)

cat("\nBias at day 90 (direction + magnitude) by scenario/n/method:\n")
b90 <- pd[pd$day == 90, c("scenario", "n", "method", "scale", "signed_bias", "rmse", "coverage")]
print(b90[order(b90$scenario, b90$n, b90$method), ], row.names = FALSE)

cat("\nM4 LOOIC form selection counts per cell:\n")
print(res$m4_selection, row.names = FALSE)

cat("\nStage 4 smoke validation complete. Outputs in ./outputs/sweep_smoke\n")
