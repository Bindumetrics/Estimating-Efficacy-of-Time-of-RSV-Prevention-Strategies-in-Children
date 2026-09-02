# ---------------------------------------------------------------------------
# Backfill the Bayesian estimator into the completed revised simulation.
#
# This reconstructs each already-simulated dataset from its original
# worker-independent random stream, fits ONLY the Bayesian estimator, verifies
# its realised event count against the stored checkpoint, then appends the
# results.  The original checkpoint is copied to outputs/checkpoints/backups/
# before it is replaced.  The script is resumable: scenarios that already have
# Bayesian rows are skipped.
#
# Usage:
#   Rscript scripts/backfill_bayes.R [workers] [scenario[,scenario,...]]
#   Rscript scripts/backfill_bayes.R 4              # all missing scenarios
#   Rscript scripts/backfill_bayes.R 4 events_30    # one scenario
#
# IMPORTANT: This is a full R=1000 backfill, not a pilot.  Do not change the
# number of replications: otherwise the added Bayesian estimates are no longer
# paired with the existing frequentist replicates.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(future); library(future.apply); library(digest); library(survival)
})
source("R/dgp.R"); source("R/estimators.R"); source("R/bayes.R")
source("R/runner.R"); source("R/scenarios.R"); source("R/metrics.R")

args <- commandArgs(trailingOnly = TRUE)
WORKERS <- if (length(args) >= 1L) as.integer(args[[1]]) else 4L
requested <- if (length(args) >= 2L) strsplit(args[[2]], ",", fixed = TRUE)[[1]] else NULL
SEED <- 20260725L

say <- function(fmt, ...) {
  cat(format(Sys.time(), "%H:%M:%S "), sprintf(fmt, ...), "\n", sep = "")
  flush.console()
}

S <- fix_ltfu_sweep_scale(build_scenarios())
if (!is.null(requested)) {
  unknown <- setdiff(requested, names(S))
  if (length(unknown)) stop("Unknown scenario(s): ", paste(unknown, collapse = ", "))
  S <- S[intersect(names(S), requested)]
}

dir.create("outputs/checkpoints/backups", recursive = TRUE, showWarnings = FALSE)
plan(multisession, workers = WORKERS)
on.exit(plan(sequential), add = TRUE)

for (scenario_name in names(S)) {
  scen <- S[[scenario_name]]
  scen_index <- match(scenario_name, names(fix_ltfu_sweep_scale(build_scenarios())))
  ck <- file.path("outputs/checkpoints", paste0(scenario_name, ".rds"))
  if (!file.exists(ck)) stop("Missing checkpoint: ", ck)
  old <- readRDS(ck)
  if ("bayes" %in% unique(old$estimates$method)) {
    say("%-22s already contains Bayesian estimates; skipping", scenario_name)
    next
  }

  reps <- sort(unique(old$estimates$rep))
  if (!identical(reps, seq_len(1000L)))
    stop(scenario_name, ": expected the completed R=1000 replicate set")

  say("%-22s reconstructing %d paired datasets and fitting Bayesian model", scenario_name, length(reps))
  t0 <- proc.time()[["elapsed"]]
  out <- future_lapply(reps, function(i) {
    run_replicate(scen, i, estimators = "bayes", bayes_seed = SEED + 1000L * scen_index + i)
  }, future.seed = SEED + scen_index,
     future.packages = c("survival", "splines", "digest", "cmdstanr", "posterior"))

  observed_events <- vapply(out, `[[`, integer(1), "n_events")
  if (!identical(observed_events, old$n_events))
    stop(scenario_name, ": reconstructed event counts differ from checkpoint; refusing to append unpaired results")

  bayes_est <- do.call(rbind, lapply(out, `[[`, "estimates"))
  bayes_dgn <- do.call(rbind, lapply(out, function(z)
    z$diagnostics[, setdiff(names(z$diagnostics), "extra"), drop = FALSE]))
  if (nrow(bayes_est) != length(reps) * nrow(old$truth) || nrow(bayes_dgn) != length(reps))
    stop(scenario_name, ": incomplete Bayesian backfill")

  updated <- old
  updated$estimates <- rbind(old$estimates, bayes_est)
  updated$diagnostics <- rbind(old$diagnostics, bayes_dgn)
  backup <- file.path("outputs/checkpoints/backups", paste0(scenario_name, "-pre-bayes.rds"))
  if (!file.exists(backup) && !file.copy(ck, backup))
    stop("Could not create checkpoint backup: ", backup)
  temporary <- paste0(ck, ".bayes-tmp")
  saveRDS(updated, temporary)
  if (!file.copy(temporary, ck, overwrite = TRUE))
    stop("Could not update checkpoint: ", ck)
  unlink(temporary)
  say("%-22s complete in %.1f min; checkpoint backed up", scenario_name,
      (proc.time()[["elapsed"]] - t0) / 60)
}

plan(sequential)
say("re-aggregating output tables and figures input")
source("scripts/aggregate.R")
say("Bayesian backfill complete")
