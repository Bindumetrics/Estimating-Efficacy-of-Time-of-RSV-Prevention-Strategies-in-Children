# validate_m3.R
# Single-dataset validation for M3 before moving to M4.

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))
source(file.path("R", "methods", "m3.R"))

config <- make_sim_config(m4_interval_width = 15L)
calibrated <- calibrate_all_scenarios(config)
scenario_id <- "B"
scenario <- calibrated[[scenario_id]]

grid <- seq.int(config$evaluation$grid_start, scenario$spec$horizon)
truth <- make_truth_grid(scenario, scenario_id, config)
dat <- simulate_trial(n = 7000L, calibrated = scenario, loss_to_followup = 0.05, seed = config$meta$rng_seed + 4L)
fit <- fit_m3(dat, grid = grid, bandwidth_days = config$m3$bandwidth_days, kernel = config$m3$kernel)

if (!isTRUE(fit$converged)) {
  failed_days <- fit$estimates[!fit$estimates$local_converged, c("day", "local_failure_reason")]
  print(head(failed_days, 10))
  stop("M3 validation fit failed: ", fit$failure_reason)
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
                            c("day", "ve", "lower", "upper", "ve_h", "ve_r", "bias_vs_ve_h", "effective_events")]

event_summary <- aggregate(status ~ z, data = dat, FUN = sum)
names(event_summary) <- c("z", "events")
event_summary$n <- as.integer(table(dat$z)[as.character(event_summary$z)])
event_summary$arm <- ifelse(event_summary$z == 1L, "vaccine", "placebo")
event_summary <- event_summary[, c("arm", "n", "events")]

write.csv(fit$estimates, file.path("outputs", "m3_validation_estimates.csv"), row.names = FALSE)
write.csv(summary_table, file.path("outputs", "m3_validation_summary_days.csv"), row.names = FALSE)
write.csv(event_summary, file.path("outputs", "m3_validation_event_summary.csv"), row.names = FALSE)
saveRDS(list(data = dat, fit = fit, comparison = comparison), file.path("outputs", "m3_validation_fit.rds"))

plot_window <- c(0, 1)
comparison$lower_plot <- pmax(plot_window[1], comparison$lower)
comparison$upper_plot <- pmin(plot_window[2], comparison$upper)

png(file.path("outputs", "m3_validation_plot.png"), width = 1400, height = 900, res = 160)
plot(comparison$day, comparison$ve_h, type = "l", lwd = 2,
     ylim = plot_window,
     xlab = "Day", ylab = "Hazard-scale VE",
     main = "M3 validation: scenario B, n = 7000, 5% loss to follow-up")
polygon(c(comparison$day, rev(comparison$day)), c(comparison$lower_plot, rev(comparison$upper_plot)),
        col = grDevices::adjustcolor("orange", alpha.f = 0.20), border = NA)
lines(comparison$day, pmax(plot_window[1], pmin(plot_window[2], comparison$ve)), lwd = 2, col = "darkorange3")
lines(comparison$day, comparison$ve_h, lwd = 2, col = "black")
mtext("CI ribbon clipped to VE range [0, 1] for display; raw interval limits are saved in CSV/RDS outputs.",
      side = 1, line = 4, cex = 0.75)
legend("bottomleft", legend = c("true VE_h(t)", "M3 estimate", "95% pointwise CI, clipped for display"),
       lty = c(1, 1, NA), lwd = c(2, 2, NA), pch = c(NA, NA, 15),
       col = c("black", "darkorange3", grDevices::adjustcolor("orange", alpha.f = 0.20)), bty = "n")
dev.off()

cat("M3 validation complete.\n")
cat("Event summary:\n")
print(event_summary)
cat("\nBandwidth:", fit$bandwidth_days, "days; kernel:", fit$kernel, "\n")
cat("\nSummary days:\n")
print(summary_table)

