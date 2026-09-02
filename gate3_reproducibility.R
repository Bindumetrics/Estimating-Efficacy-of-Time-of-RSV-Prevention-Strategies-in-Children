# ---------------------------------------------------------------------------
# GATE 3 -- reproducibility. Every check reports a NUMBER (or a hash).
#
# Parallel RNG uses L'Ecuyer-CMRG streams so results are INDEPENDENT OF WORKER
# COUNT. This is verified by running the same scenario at 2 and 4 workers and
# asserting identical output, not assumed from the fact that a seed was set.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(future); library(future.apply); library(digest); library(survival)
})
source("R/dgp.R"); source("R/estimators.R"); source("R/bayes.R")
source("R/runner.R"); source("R/scenarios.R")

say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
PASS <- TRUE

# Bayes is excluded from the reproducibility check: CmdStan runs in a separate
# process with its own RNG, so bitwise reproducibility would require pinning the
# Stan seed per replicate. That IS done in the production runner (bayes_seed is
# derived from the replicate id), but it is checked separately below rather than
# folded into the worker-count comparison.
FREQ <- c("cox_tdc", "durham", "tian", "pooled_logistic")

run_scenario_par <- function(scen, nrep, seed, workers, estimators = FREQ) {
  plan(multisession, workers = workers)
  on.exit(plan(sequential), add = TRUE)
  out <- future_lapply(seq_len(nrep), function(i) {
    run_replicate(scen, i, estimators = estimators)
  }, future.seed = seed, future.packages = c("survival", "splines", "digest"),
  future.globals = TRUE)
  do.call(rbind, lapply(out, `[[`, "estimates"))
}

scen <- make_scenario("gate3", target_events = 200)
nrep <- 24
SEED <- 20260723

say("\n=== 3.1 Worker-count independence (L'Ecuyer-CMRG streams) ===")
say("  Same scenario, same seed, %d replicates, four frequentist methods.", nrep)
hashes <- c()
for (w in c(2, 4)) {
  t0 <- proc.time()[["elapsed"]]
  r <- run_scenario_par(scen, nrep, SEED, w)
  el <- proc.time()[["elapsed"]] - t0
  h <- digest(r, algo = "sha256")
  hashes[as.character(w)] <- h
  say("  %d workers: %.1f s   sha256(estimates) = %s", w, el, substr(h, 1, 24))
}
ok31 <- length(unique(hashes)) == 1L
PASS <- PASS && ok31
say("  [%s] output identical across worker counts", if (ok31) "PASS" else "FAIL")

say("\n=== 3.2 Same-seed re-run reproduces exactly ===")
r1 <- run_scenario_par(scen, nrep, SEED, 4)
r2 <- run_scenario_par(scen, nrep, SEED, 4)
h1 <- digest(r1, algo = "sha256"); h2 <- digest(r2, algo = "sha256")
say("  run 1 sha256 = %s", substr(h1, 1, 32))
say("  run 2 sha256 = %s", substr(h2, 1, 32))
ok32 <- identical(h1, h2)
PASS <- PASS && ok32
say("  [%s] byte-identical on re-run", if (ok32) "PASS" else "FAIL")

say("\n=== 3.3 A different seed must give DIFFERENT results ===")
say("  (guards against a seed that is silently ignored -- a passing 3.2 with a")
say("   broken RNG would look identical here too)")
r3 <- run_scenario_par(scen, nrep, SEED + 1, 4)
h3 <- digest(r3, algo = "sha256")
ok33 <- !identical(h1, h3)
say("  seed+1 sha256 = %s", substr(h3, 1, 32))
PASS <- PASS && ok33
say("  [%s] different seed changes the output", if (ok33) "PASS" else "FAIL")

say("\n=== 3.4 Pairing: all estimators see byte-identical data ===")
say("  run_replicate() hashes the dataset and re-checks after every estimator;")
say("  a mutation raises an error rather than silently breaking the pairing.")
rr <- run_replicate(scen, 1L, estimators = FREQ, keep_data = TRUE)
say("  dataset hash = %s, events = %d, methods = %d",
    rr$data_hash, rr$n_events, length(unique(rr$estimates$method)))
ok34 <- identical(digest(rr$data, algo = "xxhash64"), rr$data_hash)
PASS <- PASS && ok34
say("  [%s] dataset unchanged after the full estimator loop", if (ok34) "PASS" else "FAIL")

say("\n=== 3.5 Stan seed control ===")
d <- simulate_trial(scen$n_c, scen$n_v, scen$h0, scen$ve_fn, scen$horizon,
                    theta = 0, ltfu_prob = 0.04, accrual_days = 120, scale = scen$scale)
b1 <- est_bayes(d, scen$grid, horizon = 180, seed = 555)
b2 <- est_bayes(d, scen$grid, horizon = 180, seed = 555)
b3 <- est_bayes(d, scen$grid, horizon = 180, seed = 999)
same <- isTRUE(all.equal(b1$ve_hat, b2$ve_hat))
diff <- !isTRUE(all.equal(b1$ve_hat, b3$ve_hat))
say("  same Stan seed reproduces: %s ; different seed differs: %s", same, diff)
say("  max |difference| same-seed = %.3g ; different-seed = %.3g",
    max(abs(b1$ve_hat - b2$ve_hat)), max(abs(b1$ve_hat - b3$ve_hat)))
PASS <- PASS && same && diff
say("  [%s] Stan sampling is seed-controlled", if (same && diff) "PASS" else "FAIL")

say("\n=== 3.6 Environment record ===")
dir.create("outputs/environment", showWarnings = FALSE, recursive = TRUE)
si <- sessionInfo()
writeLines(capture.output(print(si)), "outputs/environment/sessionInfo.txt")
cv <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NA_character_)
writeLines(c(paste("cmdstan_version:", cv),
             paste("R_version:", R.version.string),
             paste("seed:", SEED),
             paste("recorded:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))),
           "outputs/environment/run_metadata.txt")
say("  R: %s", R.version.string)
say("  CmdStan: %s", cv)
say("  sessionInfo -> outputs/environment/sessionInfo.txt")

say("\n=== GATE 3 SUMMARY ===")
say("  Overall: %s", if (PASS) "PASS" else "FAIL")
