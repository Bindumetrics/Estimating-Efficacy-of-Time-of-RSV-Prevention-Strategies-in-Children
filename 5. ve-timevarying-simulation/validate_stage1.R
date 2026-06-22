# validate_stage1.R
# Run from the project root with: Rscript validate_stage1.R

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))

config <- make_sim_config(m4_interval_width = 15L)
calibrated <- calibrate_all_scenarios(config)
truth_grid <- do.call(rbind, Map(make_truth_grid, calibrated, names(calibrated), MoreArgs = list(config = config)))

calibration_table <- do.call(rbind, lapply(names(calibrated), function(id) {
  x <- calibrated[[id]]
  cbind(
    scenario = id,
    label = x$spec$label,
    family = x$truth$family,
    x$calibration,
    objective = x$optimisation$objective,
    convergence = x$optimisation$convergence,
    row.names = NULL
  )
}))

incidence_table <- do.call(rbind, lapply(names(calibrated), function(id) {
  x <- calibrated[[id]]
  h <- x$spec$horizon
  risks <- cumulative_risks(h, x$baseline, x$truth)
  data.frame(
    scenario = id,
    horizon = h,
    target_placebo_risk = x$spec$target_placebo_cumulative_incidence,
    implied_placebo_risk = risks$placebo_risk,
    implied_vaccine_risk = risks$vaccine_risk,
    implied_cumulative_ve_r = 1 - risks$vaccine_risk / risks$placebo_risk,
    baseline_lambda = x$baseline$lambda
  )
}))

m4_intervals <- do.call(rbind, lapply(names(config$scenarios), function(id) {
  grid <- m4_interval_grid(config$scenarios[[id]]$horizon, config$m4$interval_width)
  cbind(scenario = id, grid)
}))

write.csv(calibration_table, file.path("outputs", "stage1_calibration_table.csv"), row.names = FALSE)
write.csv(incidence_table, file.path("outputs", "stage1_incidence_table.csv"), row.names = FALSE)
write.csv(truth_grid, file.path("outputs", "stage1_truth_grid.csv"), row.names = FALSE)
write.csv(m4_intervals, file.path("outputs", "stage1_m4_interval_grid.csv"), row.names = FALSE)
saveRDS(list(config = config, calibrated = calibrated, truth_grid = truth_grid), file.path("outputs", "stage1_calibrated_dgm.rds"))

png(file.path("outputs", "stage1_truth_curves.png"), width = 1800, height = 1200, res = 180)
op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (id in names(calibrated)) {
  one <- truth_grid[truth_grid$scenario == id, ]
  plot(one$day, one$ve_h, type = "l", ylim = c(0, 1), lwd = 2,
       xlab = "Day", ylab = "VE", main = paste(id, calibrated[[id]]$spec$label))
  lines(one$day, one$ve_r, lwd = 2, lty = 2)
  points(calibrated[[id]]$calibration$day, calibrated[[id]]$calibration$target_cumulative_ve_r,
         pch = 19, col = "firebrick")
  legend("topright", legend = c("VE_h(t)", "VE_r(0..t)", "cumulative anchors"),
         lty = c(1, 2, NA), pch = c(NA, NA, 19), bty = "n", cex = 0.75,
         col = c("black", "black", "firebrick"))
}
par(op)
dev.off()

cat("Stage 1 validation complete.\n")
cat("Outputs written to ./outputs\n")
print(calibration_table)
print(incidence_table)
