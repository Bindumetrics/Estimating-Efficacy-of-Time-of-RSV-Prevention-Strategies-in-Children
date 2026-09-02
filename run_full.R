# ---------------------------------------------------------------------------
# Full run driver. RESUMABLE: per-scenario results are written to
# outputs/checkpoints/ as they complete, and an interrupted run picks up from
# the last completed scenario rather than starting over. A thermal shutdown or
# a crash costs one scenario, not the whole grid.
#
# Do NOT run this until the Gate 4 projection has been accepted.
#
#   Rscript scripts/run_full.R [R] [workers]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(future); library(future.apply); library(digest); library(survival)
})
source("R/dgp.R"); source("R/estimators.R"); source("R/bayes.R")
source("R/runner.R"); source("R/scenarios.R"); source("R/metrics.R")

args    <- commandArgs(trailingOnly = TRUE)
R_FULL  <- if (length(args) >= 1) as.integer(args[1]) else 1000L
WORKERS <- if (length(args) >= 2) as.integer(args[2]) else 4L
SEED    <- 20260725L

# SPLIT DESIGN (agreed at the Gate 4 sign-off).
#
# The four frequentist methods cost ~0.76 s per replicate combined; the Bayesian
# method costs ~5.5 s, i.e. ~88% of the total. Running all five on all 23
# scenarios at R = 1000 projects to ~23 hours. Running the frequentist methods
# everywhere and the Bayesian method on a focused subset brings that to ~8 hours
# without weakening any frequentist comparison.
#
# The subset is chosen to span the inferentially distinct cases, not the cheap
# ones: the reference, the null (type I error), the two hardest misspecification
# shapes, the frailty scenario where every method is misspecified, and the
# verified 360-day applied anchor.
#
# SCOPE LIMIT TO REPORT: Bayesian comparisons cover these six scenarios only.
# Any claim about the Bayesian method's behaviour under the sample-size, dropout
# or baseline-hazard sweeps is NOT supported by this run.
BAYES_SUBSET <- c("ref", "shape_constant", "shape_rapid", "shape_piecewise",
                  "frailty_1.0", "applied_matisse_360")

# Pairing is preserved across the split: each scenario runs ONE future_lapply
# over replicates, and the estimator list is chosen inside it. Replicate i of a
# scenario therefore generates one dataset that every method run on it shares,
# exactly as in the non-split design.
methods_for <- function(scenario_name) {
  if (scenario_name %in% BAYES_SUBSET) ESTIMATORS else setdiff(ESTIMATORS, "bayes")
}

say <- function(fmt, ...) {
  cat(format(Sys.time(), "%H:%M:%S "), sprintf(fmt, ...), "\n", sep = "")
  flush.console()
}

dir.create("outputs/checkpoints", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

S <- fix_ltfu_sweep_scale(build_scenarios())
write.csv(scenario_register(S), "outputs/scenario_register.csv", row.names = FALSE)
say("%d scenarios, R = %d, %d workers, seed = %d", length(S), R_FULL, WORKERS, SEED)

plan(multisession, workers = WORKERS)
on.exit(plan(sequential), add = TRUE)

t_start <- proc.time()[["elapsed"]]
for (k in seq_along(S)) {
  scen <- S[[k]]
  ck <- file.path("outputs/checkpoints", paste0(scen$name, ".rds"))
  if (file.exists(ck)) { say("[%2d/%2d] %-22s already done, skipping", k, length(S), scen$name); next }

  t0 <- proc.time()[["elapsed"]]
  meth <- methods_for(scen$name)
  # future.seed = a single integer gives reproducible L'Ecuyer-CMRG streams that
  # are independent of the worker count (verified at Gate 3.1). Offsetting by the
  # scenario index keeps scenarios from sharing a stream while remaining
  # deterministic.
  out <- future_lapply(seq_len(R_FULL), function(i) {
    run_replicate(scen, i, estimators = meth, bayes_seed = SEED + 1000L * k + i)
  }, future.seed = SEED + k,
     future.packages = c("survival", "splines", "digest", "cmdstanr", "posterior"))

  est <- do.call(rbind, lapply(out, `[[`, "estimates"))
  dgn <- do.call(rbind, lapply(out, function(z)
    z$diagnostics[, setdiff(names(z$diagnostics), "extra")]))
  truth <- data.frame(t = scen$grid, ve_true = scen$ve_h_true,
                      ve_id_true = scen$ve_id_true)

  saveRDS(list(scenario = scen$name, estimates = est, diagnostics = dgn,
               truth = truth, n_events = sapply(out, `[[`, "n_events"),
               scale = scen$scale, seed = SEED + k), ck)

  el <- proc.time()[["elapsed"]] - t0
  fail <- tapply(!dgn$converged, dgn$method, mean)
  say("[%2d/%2d] %-22s %6.1f s (%.2f s/rep)  events %.0f  %s  failures: %s",
      k, length(S), scen$name, el, el / R_FULL,
      median(sapply(out, `[[`, "n_events")),
      if ("bayes" %in% meth) "[+bayes]" else "[freq  ]",
      paste(sprintf("%s %.0f%%", names(fail), 100 * fail), collapse = " "))
}
plan(sequential)
say("all scenarios complete in %.2f hours", (proc.time()[["elapsed"]] - t_start) / 3600)

# ---------------------------------------------------------------------------
say("aggregating")
cks <- list.files("outputs/checkpoints", pattern = "\\.rds$", full.names = TRUE)
cks <- cks[basename(cks) != "pilot.rds"]
all_est <- do.call(rbind, lapply(cks, function(f) readRDS(f)$estimates))
all_dgn <- do.call(rbind, lapply(cks, function(f) readRDS(f)$diagnostics))
all_tru <- do.call(rbind, lapply(cks, function(f) {
  z <- readRDS(f); cbind(scenario = z$scenario, z$truth) }))

# Raw estimates are persisted so every table regenerates without re-running.
if (requireNamespace("arrow", quietly = TRUE)) {
  arrow::write_parquet(all_est, "outputs/raw_estimates.parquet")
} else {
  saveRDS(all_est, "outputs/raw_estimates.rds")
  say("NOTE: arrow unavailable; raw estimates written as .rds instead of .parquet")
}
write.csv(all_dgn, "outputs/replicate_diagnostics.csv", row.names = FALSE)

pw <- do.call(rbind, lapply(split(all_est, all_est$scenario), function(g) {
  tr <- all_tru[all_tru$scenario == g$scenario[1], c("t", "ve_true")]
  pointwise_metrics(g, tr)
}))
write.csv(pw, "outputs/pointwise_metrics.csv", row.names = FALSE)

sm <- summary_metrics(pw)
write.csv(sm, "outputs/summary_metrics.csv", row.names = FALSE)

pc <- do.call(rbind, lapply(split(all_est, all_est$scenario), function(g) {
  tr <- all_tru[all_tru$scenario == g$scenario[1], c("t", "ve_true")]
  paired_comparisons(g, tr)
}))
write.csv(pc, "outputs/paired_comparisons.csv", row.names = FALSE)

pt <- power_typeI(all_dgn)
if (!is.null(pt)) write.csv(pt, "outputs/power_typeI.csv", row.names = FALSE)

rk <- rank_matrix(sm)
write.csv(rk, "outputs/rank_matrix.csv", row.names = FALSE)
say("rank matrix reproduces from its printed columns: %s",
    attr(rk, "mean_rank_reproduces"))

bs <- do.call(rbind, lapply(split(all_est, all_est$scenario), function(g) {
  tr <- all_tru[all_tru$scenario == g$scenario[1], c("t", "ve_true")]
  blowup_sensitivity(g, tr)
}))
write.csv(bs, "outputs/blowup_sensitivity.csv", row.names = FALSE)

st <- curve_stability(all_est)
write.csv(aggregate(mean_abs_d2 ~ scenario + method, st, mean),
          "outputs/curve_stability.csv", row.names = FALSE)

say("done. %d rows of raw estimates across %d scenarios.",
    nrow(all_est), length(unique(all_est$scenario)))
say("Inconclusive paired comparisons (|diff| < 1.96 MCSE): %d of %d",
    sum(!pc$conclusive), nrow(pc))
