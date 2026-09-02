# ---------------------------------------------------------------------------
# Aggregate checkpointed results into the output tables. Reads only the
# per-scenario .rds files, so it is cheap and can be re-run without touching the
# simulation. Separated from run_full.R precisely so an aggregation bug never
# costs the run again.
#
#   Rscript scripts/aggregate.R
# ---------------------------------------------------------------------------

source("R/dgp.R"); source("R/metrics.R")
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

cks <- list.files("outputs/checkpoints", pattern = "\\.rds$", full.names = TRUE)
cks <- cks[basename(cks) != "pilot.rds"]
say("aggregating %d scenario checkpoints", length(cks))

all_est <- do.call(rbind, lapply(cks, function(f) readRDS(f)$estimates))
all_dgn <- do.call(rbind, lapply(cks, function(f) readRDS(f)$diagnostics))
all_tru <- do.call(rbind, lapply(cks, function(f) {
  z <- readRDS(f); cbind(scenario = z$scenario, z$truth) }))

if (requireNamespace("arrow", quietly = TRUE)) {
  arrow::write_parquet(all_est, "outputs/raw_estimates.parquet")
  say("wrote outputs/raw_estimates.parquet (%d rows)", nrow(all_est))
} else {
  saveRDS(all_est, "outputs/raw_estimates.rds")
}
write.csv(all_dgn, "outputs/replicate_diagnostics.csv", row.names = FALSE)

pw <- do.call(rbind, lapply(split(all_est, all_est$scenario), function(g) {
  tr <- unique(all_tru[all_tru$scenario == g$scenario[1], c("t", "ve_true")])
  pointwise_metrics(g, tr)
}))
write.csv(pw, "outputs/pointwise_metrics.csv", row.names = FALSE)
say("pointwise_metrics: %d rows", nrow(pw))

sm <- summary_metrics(pw[!is.na(pw$bias), ])
write.csv(sm, "outputs/summary_metrics.csv", row.names = FALSE)

pc <- do.call(rbind, lapply(split(all_est, all_est$scenario), function(g) {
  tr <- unique(all_tru[all_tru$scenario == g$scenario[1], c("t", "ve_true")])
  paired_comparisons(g, tr)
}))
write.csv(pc, "outputs/paired_comparisons.csv", row.names = FALSE)
say("paired_comparisons: %d rows, %d inconclusive (|diff| < 1.96 MCSE)",
    nrow(pc), sum(!pc$conclusive))

pt <- power_typeI(all_dgn)
if (!is.null(pt)) { write.csv(pt, "outputs/power_typeI.csv", row.names = FALSE)
  say("power_typeI: %d rows", nrow(pt)) }

rk <- rank_matrix(sm)
write.csv(rk, "outputs/rank_matrix.csv", row.names = FALSE)
say("rank matrix reproduces from its printed columns: %s",
    attr(rk, "mean_rank_reproduces"))

bs <- do.call(rbind, lapply(split(all_est, all_est$scenario), function(g) {
  tr <- unique(all_tru[all_tru$scenario == g$scenario[1], c("t", "ve_true")])
  blowup_sensitivity(g, tr)
}))
write.csv(bs, "outputs/blowup_sensitivity.csv", row.names = FALSE)

st <- curve_stability(all_est)
write.csv(aggregate(mean_abs_d2 ~ scenario + method, st, mean),
          "outputs/curve_stability.csv", row.names = FALSE)

say("done. %d rows of raw estimates across %d scenarios, %d methods.",
    nrow(all_est), length(unique(all_est$scenario)), length(unique(all_est$method)))
