# validate_m4.R
# Single-dataset validation for M4 (Bayesian binomial interval-count model).
# Fits all three candidate forms, reports LOOIC selection, and compares the
# selected posterior VE(t) curve to the DGM truth.
#
# Run from project root:  Rscript validate_m4.R

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))
source(file.path("R", "methods", "m4.R"))

config <- make_sim_config(m4_interval_width = 15L)
calibrated <- calibrate_all_scenarios(config)
scenario_id <- "B"
scenario <- calibrated[[scenario_id]]
horizon <- scenario$spec$horizon

grid <- seq.int(config$evaluation$grid_start, horizon)
truth <- make_truth_grid(scenario, scenario_id, config)

# DGM-implied interval-conditional risk-scale VE: the quantity M4 actually
# targets (1 - conditional vaccine risk / conditional placebo risk per interval),
# mapped onto the day grid. This is the fair target for M4, distinct from the
# cumulative VE_r and from the hazard-scale VE_h.
interval_grid <- m4_interval_grid(horizon, config$m4$interval_width)
ve_r_interval_truth <- function(days) {
  H_c <- function(t) baseline_cumhaz(t, scenario$baseline)
  H_v <- function(t) vaccine_cumhaz(t, scenario$baseline, scenario$truth)
  val <- numeric(length(days))
  for (i in seq_along(days)) {
    d <- days[i]
    k <- which(interval_grid$start < d & interval_grid$end >= d)
    if (length(k) == 0) k <- nrow(interval_grid)
    a <- interval_grid$start[k]; b <- interval_grid$end[k]
    pi_p <- 1 - exp(-(H_c(b) - H_c(a)))
    pi_v <- 1 - exp(-(H_v(b) - H_v(a)))
    val[i] <- 1 - pi_v / pi_p
  }
  val
}
truth$ve_r_interval <- ve_r_interval_truth(truth$day)

dat <- simulate_trial(n = 7000L, calibrated = scenario, loss_to_followup = 0.05,
                      seed = config$meta$rng_seed + 5L)

model <- m4_stan_model()
fit <- fit_m4(dat, grid = grid, horizon = horizon,
              interval_width = config$m4$interval_width, model = model,
              seed = config$meta$rng_seed + 5L,
              iter_warmup = 1000L, iter_sampling = 1000L, parallel_chains = 4L)

if (!isTRUE(fit$converged)) stop("M4 validation fit failed: ", fit$failure_reason)

cat("LOOIC table (lower is better):\n")
print(fit$looic_table)
cat("\nLOOIC-selected form:", fit$selected_form, "\n\n")

comparison <- merge(fit$estimates, truth[, c("day", "ve_h", "ve_r", "ve_r_interval")],
                    by = "day", all.x = TRUE)
comparison$bias_vs_ve_r_interval <- comparison$ve - comparison$ve_r_interval

summary_days <- config$evaluation$summary_days[config$evaluation$summary_days <= horizon]
summary_table <- comparison[comparison$day %in% summary_days,
                            c("day", "ve", "lower", "upper", "ve_h", "ve_r",
                              "ve_r_interval", "bias_vs_ve_r_interval")]

event_summary <- aggregate(status ~ z, data = dat, FUN = sum)
names(event_summary) <- c("z", "events")
event_summary$n <- as.integer(table(dat$z)[as.character(event_summary$z)])
event_summary$arm <- ifelse(event_summary$z == 1L, "vaccine", "placebo")
event_summary <- event_summary[, c("arm", "n", "events")]

write.csv(fit$looic_table, file.path("outputs", "m4_validation_looic.csv"), row.names = FALSE)
write.csv(fit$all_form_estimates, file.path("outputs", "m4_validation_all_forms.csv"), row.names = FALSE)
write.csv(summary_table, file.path("outputs", "m4_validation_summary_days.csv"), row.names = FALSE)
write.csv(fit$bin, file.path("outputs", "m4_validation_bin_counts.csv"), row.names = FALSE)
write.csv(event_summary, file.path("outputs", "m4_validation_event_summary.csv"), row.names = FALSE)
saveRDS(list(data = dat, fit = fit, comparison = comparison), file.path("outputs", "m4_validation_fit.rds"))

png(file.path("outputs", "m4_validation_plot.png"), width = 1400, height = 900, res = 160)
plot(comparison$day, comparison$ve_r_interval, type = "l", lwd = 2,
     ylim = c(0, 1), xlab = "Day", ylab = "Risk-scale VE",
     main = sprintf("M4 validation: scenario B, n = 7000 (LOOIC form: %s)", fit$selected_form))
polygon(c(comparison$day, rev(comparison$day)), c(comparison$lower, rev(comparison$upper)),
        col = grDevices::adjustcolor("seagreen", alpha.f = 0.20), border = NA)
lines(comparison$day, comparison$ve, lwd = 2, col = "seagreen4")
lines(comparison$day, comparison$ve_r_interval, lwd = 2, col = "black")
lines(comparison$day, comparison$ve_h, lwd = 1, lty = 3, col = "grey40")
legend("bottomleft",
       legend = c("true interval VE_r(t) (M4 target)", "M4 selected estimate",
                  "95% CrI", "true VE_h(t) (hazard scale)"),
       lty = c(1, 1, NA, 3), lwd = c(2, 2, NA, 1), pch = c(NA, NA, 15, NA),
       col = c("black", "seagreen4", grDevices::adjustcolor("seagreen", alpha.f = 0.20), "grey40"),
       bty = "n")
dev.off()

cat("M4 validation complete.\n")
cat("Event summary:\n"); print(event_summary)
cat("\nInterval width:", fit$interval_width, "days;  intervals:", nrow(fit$bin), "\n")
cat("\nSummary days:\n"); print(summary_table)
