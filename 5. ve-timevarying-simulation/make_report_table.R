# make_report_table.R
# Collapse the quick-run summary tables into a single per-method reporting table
# for M1-M3 (M1b shown as supplementary). All on the hazard scale, so directly
# comparable. Averages over scenarios x n x loss x reporting-days.

sm <- file.path("outputs", "summary")
day_cols <- c("d30", "d60", "d90", "d120", "d150", "d180")

bias <- read.csv(file.path(sm, "bias.csv"))      # signed_bias by day
rmse <- read.csv(file.path(sm, "rmse.csv"))
cov  <- read.csv(file.path(sm, "coverage.csv"))
wid  <- read.csv(file.path(sm, "width.csv"))
ise  <- read.csv(file.path(sm, "ise.csv"))
reff <- read.csv(file.path(sm, "relative_efficiency.csv"))  # hazard-scale ranks

# Keep the hazard-scale rows only (M1-M3, M1b live on the hazard scale).
haz <- function(d) if ("scale" %in% names(d)) d[d$scale == "hazard", ] else d

row_mean <- function(d) rowMeans(d[, day_cols], na.rm = TRUE)

# Per (cell, method) day-averaged values, then MEDIAN across cells per method.
# Medians are used because a handful of cells produce numerically degenerate
# (blown-up) estimates at R=10; means are dominated by those outliers.
agg_metric <- function(d, fn = median, absval = FALSE) {
  d <- haz(d)
  v <- row_mean(d)
  if (absval) v <- abs(v)
  tapply(v, d$method, fn, na.rm = TRUE)
}

med_abs_bias <- agg_metric(bias, absval = TRUE)
med_signed   <- agg_metric(bias, absval = FALSE)
med_rmse     <- agg_metric(rmse)
mean_cov     <- agg_metric(cov, fn = mean)   # coverage is bounded [0,1]; mean is fine
med_width    <- agg_metric(wid)

ise_h <- haz(ise)
med_ise <- tapply(ise_h$mean_ise, ise_h$method, median, na.rm = TRUE)

mean_rank <- tapply(reff$rank, reff$method, mean, na.rm = TRUE)  # rank is bounded 1-4

# Instability: fraction of cell-day bias values that are a gross failure
# (|bias| > 1 on the VE scale, i.e. the estimate has diverged).
bias_h <- haz(bias)
blow <- tapply(seq_len(nrow(bias_h)), bias_h$method, function(idx) {
  vals <- as.matrix(bias_h[idx, day_cols])
  mean(abs(vals) > 1, na.rm = TRUE)
})

methods <- c("M1", "M2", "M3", "M1b_supplementary")
methods <- methods[methods %in% names(med_abs_bias)]

tab <- data.frame(
  method          = methods,
  median_abs_bias = round(med_abs_bias[methods], 4),
  median_signed_bias = round(med_signed[methods], 4),
  median_rmse     = round(med_rmse[methods], 4),
  coverage_95     = round(mean_cov[methods], 3),
  median_ci_width = round(med_width[methods], 4),
  median_ise      = round(med_ise[methods], 5),
  mean_rank_1to4  = round(mean_rank[methods], 2),
  blowup_rate     = round(blow[methods], 3),
  row.names = NULL
)

cat("\n==== M1-M3 REPORTING TABLE (hazard-scale VE_h; medians over all cells & days) ====\n")
cat("Source: quick run, 5 scenarios x {500,2000,5000} x 4 loss, R=10\n")
cat("blowup_rate = fraction of cell-days with |bias| > 1 (diverged estimate)\n\n")
print(tab, row.names = FALSE)

write.csv(tab, file.path(sm, "method_report_table.csv"), row.names = FALSE)
cat("\nSaved:", file.path(sm, "method_report_table.csv"), "\n")

# Per-scenario best (lowest MEDIAN RMSE among M1-M3) for context.
cat("\n==== Best method per scenario (lowest median RMSE, M1-M3 core) ====\n")
core <- haz(rmse)
core <- core[core$method %in% c("M1", "M2", "M3"), ]
core$rmse_avg <- row_mean(core)
by_sc <- aggregate(rmse_avg ~ scenario + method, data = core, FUN = median)
for (s in sort(unique(by_sc$scenario))) {
  d <- by_sc[by_sc$scenario == s, ]
  best <- d$method[which.min(d$rmse_avg)]
  cat(sprintf("  Scenario %s: %s  (RMSE %.4f)\n", s, best, min(d$rmse_avg)))
}
