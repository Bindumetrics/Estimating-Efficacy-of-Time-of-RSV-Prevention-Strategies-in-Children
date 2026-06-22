# validate_m1.R
# Single-dataset validation for M1 before moving to the next method.

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))
source(file.path("R", "methods", "m1.R"))

config <- make_sim_config(m4_interval_width = 15L)
calibrated <- calibrate_all_scenarios(config)
scenario_id <- "B"
scenario <- calibrated[[scenario_id]]

grid <- seq.int(config$evaluation$grid_start, scenario$spec$horizon)
truth <- make_truth_grid(scenario, scenario_id, config)
set.seed(config$meta$rng_seed)
dat <- simulate_trial(n = 7000L, calibrated = scenario, loss_to_followup = 0.05, seed = config$meta$rng_seed + 1L)
fit <- fit_m1(dat, grid = grid)

if (!isTRUE(fit$converged)) {
  stop("M1 validation fit failed: ", fit$failure_reason)
}

comparison <- merge(
  fit$estimates,
  truth[, c("day", "ve_h", "ve_r")],
  by = "day",
  all.x = TRUE
)
comparison$bias_vs_ve_h <- comparison$ve - comparison$ve_h

summary_days <- config$evaluation$summary_days[config$evaluation$summary_days <= scenario$spec$horizon]
summary_table <- comparison[comparison$day %in% summary_days,
                            c("day", "ve", "lower", "upper", "ve_h", "ve_r", "bias_vs_ve_h")]

event_summary <- aggregate(status ~ z, data = dat, FUN = sum)
names(event_summary) <- c("z", "events")
event_summary$n <- as.integer(table(dat$z)[as.character(event_summary$z)])
event_summary$arm <- ifelse(event_summary$z == 1L, "vaccine", "placebo")
event_summary <- event_summary[, c("arm", "n", "events")]

write.csv(fit$estimates, file.path("outputs", "m1_validation_estimates.csv"), row.names = FALSE)
write.csv(summary_table, file.path("outputs", "m1_validation_summary_days.csv"), row.names = FALSE)
write.csv(event_summary, file.path("outputs", "m1_validation_event_summary.csv"), row.names = FALSE)
saveRDS(list(data = dat, fit = fit, comparison = comparison), file.path("outputs", "m1_validation_fit.rds"))

png(file.path("outputs", "m1_validation_plot.png"), width = 1400, height = 900, res = 160)
plot(comparison$day, comparison$ve_h, type = "l", lwd = 2, ylim = range(c(comparison$lower, comparison$upper, comparison$ve_h), finite = TRUE),
     xlab = "Day", ylab = "Hazard-scale VE", main = "M1 validation: scenario B, n = 7000, 5% loss to follow-up")
polygon(c(comparison$day, rev(comparison$day)), c(comparison$lower, rev(comparison$upper)),
        col = grDevices::adjustcolor("steelblue", alpha.f = 0.20), border = NA)
lines(comparison$day, comparison$ve, lwd = 2, col = "steelblue")
lines(comparison$day, comparison$ve_h, lwd = 2, col = "black")
legend("topright", legend = c("true VE_h(t)", "M1 estimate", "95% pointwise CI"),
       lty = c(1, 1, NA), lwd = c(2, 2, NA), pch = c(NA, NA, 15),
       col = c("black", "steelblue", grDevices::adjustcolor("steelblue", alpha.f = 0.20)), bty = "n")
dev.off()

cat("M1 validation complete.\n")
cat("Event summary:\n")
print(event_summary)
cat("\nSummary days:\n")
print(summary_table)
