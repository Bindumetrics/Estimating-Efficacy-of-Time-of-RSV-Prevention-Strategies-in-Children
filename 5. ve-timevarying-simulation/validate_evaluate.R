# validate_evaluate.R
# Stage 3 validation: fit all methods on one simulated dataset, harmonise onto
# the common day grid, score each against the correct truth scale, and report
# the full metric set. This exercises the evaluation machinery end-to-end before
# the sweep; per-cell metrics here are from a single replication.
#
# Run from project root:  Rscript validate_evaluate.R

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))
source(file.path("R", "evaluate.R"))
source(file.path("R", "methods", "m1.R"))
source(file.path("R", "methods", "m1b.R"))
source(file.path("R", "methods", "m2.R"))
source(file.path("R", "methods", "m3.R"))
source(file.path("R", "methods", "m4.R"))

config <- make_sim_config(m4_interval_width = 15L)
calibrated <- calibrate_all_scenarios(config)
scenario_id <- "B"
scenario <- calibrated[[scenario_id]]
horizon <- scenario$spec$horizon
grid <- seq.int(config$evaluation$grid_start, horizon)

n <- 5000L
ltfu <- 0.05
dat <- simulate_trial(n = n, calibrated = scenario, loss_to_followup = ltfu,
                      seed = config$meta$rng_seed + 6L)

cat("Fitting methods on one dataset (scenario B, n =", n, ")...\n")
model <- m4_stan_model()
fits <- list(
  m1  = fit_m1(dat, grid = grid),
  m1b = fit_m1b(dat, grid = grid),
  m2  = fit_m2(dat, grid = grid, basis_dimension = config$m2$basis_dimension),
  m3  = fit_m3(dat, grid = grid, bandwidth_days = config$m3$bandwidth_days, kernel = config$m3$kernel),
  m4  = fit_m4(dat, grid = grid, horizon = horizon, interval_width = config$m4$interval_width,
               model = model, seed = config$meta$rng_seed + 6L, refresh = 0L)
)
for (nm in names(fits)) {
  f <- fits[[nm]]
  cat(sprintf("  %-4s converged=%s  %s\n", f$method, isTRUE(f$converged),
              if (isTRUE(f$converged)) "" else paste0("(", f$failure_reason, ")")))
}

ev <- evaluate_methods(fits, scenario, scenario_id, config,
                       replication = 1L, n = n, loss_to_followup = ltfu)

per_day <- aggregate_per_day(ev$scored)
convergence <- aggregate_convergence(ev$scored)
ise_tab <- aggregate_ise(ev$ise)
re_hazard <- relative_efficiency(per_day, scenario_id, n, ltfu, scale = "hazard")

summary_days <- config$evaluation$summary_days[config$evaluation$summary_days <= horizon]
summary_table <- per_day[per_day$day %in% summary_days,
                         c("method", "scale", "day", "truth", "signed_bias",
                           "rmse", "coverage", "mean_width")]
summary_table <- summary_table[order(summary_table$method, summary_table$day), ]

write.csv(ev$scored, file.path("outputs", "eval_validation_per_day.csv"), row.names = FALSE)
write.csv(summary_table, file.path("outputs", "eval_validation_summary_days.csv"), row.names = FALSE)
write.csv(ise_tab, file.path("outputs", "eval_validation_ise.csv"), row.names = FALSE)
write.csv(convergence, file.path("outputs", "eval_validation_convergence.csv"), row.names = FALSE)
if (!is.null(re_hazard)) {
  write.csv(as.data.frame(re_hazard$ratio), file.path("outputs", "eval_validation_relative_efficiency.csv"))
}
saveRDS(list(fits = fits, evaluation = ev, per_day = per_day), file.path("outputs", "eval_validation.rds"))

# Combined plot: both truth scales + each method's estimate on the day grid.
targets <- truth_targets(scenario, scenario_id, config)
method_colours <- c(M1 = "steelblue", M1b_supplementary = "purple", M2 = "darkorange3",
                    M3 = "firebrick", M4 = "seagreen4")
png(file.path("outputs", "eval_validation_plot.png"), width = 1500, height = 950, res = 160)
plot(targets$day, targets$truth_hazard, type = "l", lwd = 3, col = "black", ylim = c(0, 1),
     xlab = "Day", ylab = "VE(t)",
     main = sprintf("Stage 3 validation: all methods vs truth (scenario B, n = %d)", n))
lines(targets$day, targets$truth_risk, lwd = 3, col = "black", lty = 2)
for (m in names(method_colours)) {
  s <- ev$scored[ev$scored$method == m, ]
  if (nrow(s) > 0) lines(s$day, pmin(1, pmax(-0.2, s$ve)), lwd = 2, col = method_colours[[m]])
}
legend("bottomleft", bty = "n", cex = 0.8,
       legend = c("truth VE_h (M1-M3 target)", "truth interval VE_r (M4 target)",
                  "M1", "M1b (supplementary)", "M2", "M3", "M4"),
       lwd = c(3, 3, 2, 2, 2, 2, 2),
       lty = c(1, 2, 1, 1, 1, 1, 1),
       col = c("black", "black", method_colours[c("M1", "M1b_supplementary", "M2", "M3", "M4")]))
dev.off()

cat("\nSummary-day metrics (single replication):\n")
print(summary_table, row.names = FALSE)
cat("\nISE (lower is better; M4 on risk scale, others on hazard scale):\n")
print(ise_tab[, c("method", "scale", "mean_ise")], row.names = FALSE)
cat("\nConvergence:\n")
print(convergence[, c("method", "convergence_rate", "mean_frac_days_available")], row.names = FALSE)
if (!is.null(re_hazard)) {
  cat("\nRelative efficiency (integrated MSE ratio, hazard-scale methods; row / column):\n")
  print(round(re_hazard$ratio, 3))
}
cat("\nStage 3 validation complete. Outputs written to ./outputs/eval_validation_*\n")
