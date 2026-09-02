# ---------------------------------------------------------------------------
# GATE 4 -- pilot, profiling, and a runtime projection for the full run.
#
# The brief (sec 2, sec 5) requires: benchmark worker counts rather than assume
# one; profile where time actually goes; estimate the Monte Carlo SE of the
# primary metric and CHOOSE R from it; and report a projected wall-clock time
# BEFORE launching the full run.
#
# Run: Rscript scripts/gate4_pilot.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(future); library(future.apply); library(digest); library(survival)
})
source("R/dgp.R"); source("R/estimators.R"); source("R/bayes.R")
source("R/runner.R"); source("R/scenarios.R")

say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
dir.create("outputs/environment", showWarnings = FALSE, recursive = TRUE)

S <- fix_ltfu_sweep_scale(build_scenarios())
say("\n=== 4.0 Scenario register ===")
reg <- scenario_register(S)
write.csv(reg, "outputs/scenario_register.csv", row.names = FALSE)
say("  %d scenarios (%d verified, %d assumed) -> outputs/scenario_register.csv",
    nrow(reg), sum(reg$tag == "verified"), sum(reg$tag == "assumed"))
print(reg[, c("scenario", "n_c", "n_v", "horizon", "target_events",
              "ve_shape", "theta", "ltfu_prob", "tag")], row.names = FALSE)

# ---------------------------------------------------------------------------
say("\n=== 4.1 Per-method cost profile (single replicate, reference scenario) ===")
scen <- S[["ref"]]
NPROF <- 12
tm <- matrix(NA, NPROF, length(ESTIMATORS), dimnames = list(NULL, ESTIMATORS))
for (i in seq_len(NPROF)) {
  rr <- run_replicate(scen, i, bayes_seed = 1000 + i)
  tm[i, ] <- rr$diagnostics$time_sec[match(ESTIMATORS, rr$diagnostics$method)]
}
say("  method            median s     IQR s      share of total")
tot <- sum(apply(tm, 2, median, na.rm = TRUE))
for (m in ESTIMATORS) {
  q <- quantile(tm[, m], c(0.25, 0.5, 0.75), na.rm = TRUE)
  say("  %-16s %8.3f   %6.3f-%.3f   %5.1f%%", m, q[2], q[1], q[3],
      100 * q[2] / tot)
}
say("  total per replicate (median, serial) = %.2f s", tot)
write.csv(as.data.frame(tm), "outputs/environment/timing_profile.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
say("\n=== 4.2 Worker-count benchmark ===")
say("  4 physical / 8 logical cores. Hyperthreads contend for cache on")
say("  compute-bound work, so the optimum is MEASURED, not assumed to be 7.")
bench_n <- 16
for (w in c(1, 3, 4, 6)) {
  plan(multisession, workers = w)
  t0 <- proc.time()[["elapsed"]]
  invisible(future_lapply(seq_len(bench_n), function(i) run_replicate(scen, i, bayes_seed = 5000 + i),
                          future.seed = 42,
                          future.packages = c("survival", "splines", "digest",
                                              "cmdstanr", "posterior")))
  el <- proc.time()[["elapsed"]] - t0
  say("  %d workers: %6.1f s for %d replicates  (%.2f s/rep, speedup %.2fx)",
      w, el, bench_n, el / bench_n, (tot * bench_n) / el)
  if (w == 1) base_el <- el
  assign(paste0("el_", w), el)
}
plan(sequential)

# ---------------------------------------------------------------------------
say("\n=== 4.3 Sustained-throughput check (thermal throttling) ===")
say("  U-series parts throttle: a 60-second benchmark can overstate sustained")
say("  throughput by 20-30%%. Measured over >3 minutes at the chosen worker count.")
WORKERS <- 4L
plan(multisession, workers = WORKERS)
t0 <- proc.time()[["elapsed"]]; done <- 0; marks <- c()
while (proc.time()[["elapsed"]] - t0 < 200) {
  invisible(future_lapply(1:WORKERS, function(i) run_replicate(scen, done + i, bayes_seed = 7000 + done + i),
                          future.seed = 99,
                          future.packages = c("survival", "splines", "digest",
                                              "cmdstanr", "posterior")))
  done <- done + WORKERS
  marks <- c(marks, (proc.time()[["elapsed"]] - t0) / done)
}
plan(sequential)
say("  %d replicates in %.0f s", done, proc.time()[["elapsed"]] - t0)
say("  s/replicate: first quarter %.3f, last quarter %.3f  (drift %+.1f%%)",
    mean(head(marks, max(1, length(marks) %/% 4))),
    mean(tail(marks, max(1, length(marks) %/% 4))),
    100 * (mean(tail(marks, max(1, length(marks) %/% 4))) /
             mean(head(marks, max(1, length(marks) %/% 4))) - 1))
sustained <- tail(marks, 1)
say("  sustained rate used for projection: %.3f s/replicate", sustained)

# ---------------------------------------------------------------------------
say("\n=== 4.4 Pilot: 50 replicates on three scenarios ===")
pilot_names <- c("ref", "shape_constant", "frailty_1.0")
plan(multisession, workers = WORKERS)
pilot <- list()
for (nm in pilot_names) {
  sc <- S[[nm]]
  out <- future_lapply(1:50, function(i) run_replicate(sc, i, bayes_seed = 20000 + i),
                       future.seed = 4242,
                       future.packages = c("survival", "splines", "digest",
                                           "cmdstanr", "posterior"))
  est <- do.call(rbind, lapply(out, `[[`, "estimates"))
  dgn <- do.call(rbind, lapply(out, function(z)
    z$diagnostics[, setdiff(names(z$diagnostics), "extra")]))
  pilot[[nm]] <- list(est = est, dgn = dgn, truth = sc$ve_h_true, grid = sc$grid)
  say("\n  -- %s --", nm)
  say("     events per replicate: median %.0f (IQR %.0f-%.0f)",
      median(sapply(out, `[[`, "n_events")),
      quantile(sapply(out, `[[`, "n_events"), 0.25),
      quantile(sapply(out, `[[`, "n_events"), 0.75))
  for (m in ESTIMATORS) {
    dd <- dgn[dgn$method == m, ]
    ee <- est[est$method == m, ]
    say("     %-16s converged %3.0f%%  NA cells %4.1f%%  median %.3f s  notes: %s",
        m, 100 * mean(dd$converged), 100 * mean(is.na(ee$ve_hat)),
        median(dd$time_sec, na.rm = TRUE),
        { u <- unique(dd$note[nzchar(dd$note)]); if (!length(u)) "-" else
          paste0(substr(u[1], 1, 44), if (length(u) > 1) sprintf(" (+%d more)", length(u) - 1) else "") })
  }
}
plan(sequential)
saveRDS(pilot, "outputs/checkpoints/pilot.rds")

# ---------------------------------------------------------------------------
say("\n=== 4.5 Choosing R from the pilot MCSE ===")
say("  Coverage is the binding constraint: MCSE = sqrt(p(1-p)/R).")
say("     R      MCSE(coverage) at p=0.95   half-width of 95%% MC interval")
for (R in c(250, 500, 1000, 2000)) {
  s <- sqrt(0.95 * 0.05 / R)
  say("  %5d          %5.2f%%                      +/- %.2f%%", R, 100 * s, 100 * 1.96 * s)
}
pe <- pilot[["ref"]]$est
emp_sd <- sapply(ESTIMATORS, function(m) {
  x <- pe[pe$method == m, ]
  median(tapply(x$ve_hat, x$t, sd, na.rm = TRUE), na.rm = TRUE)
})
say("\n  Pilot empirical SD of VE_hat (median over grid), by method:")
for (m in ESTIMATORS) say("    %-16s %.4f", m, emp_sd[m])
say("\n  MCSE(bias) = empSD/sqrt(R). To resolve a between-method bias gap of")
say("  0.01 VE units at 2 MCSE, R must satisfy R >= (2*empSD/0.01)^2:")
for (m in ESTIMATORS) say("    %-16s R >= %6.0f", m, (2 * emp_sd[m] / 0.01)^2)
say("\n  RECOMMENDATION: R = 1000. Coverage MCSE 0.69%% resolves a 3-point")
say("  coverage gap; going to R = 1900 for 0.50%% roughly doubles a CPU-bound")
say("  run on a 4-core U-series laptop for little inferential gain.")

# ---------------------------------------------------------------------------
say("\n=== 4.6 Full-run projection ===")
R_FULL <- 1000
n_scen <- length(S)
secs <- sustained * R_FULL * n_scen
say("  %d scenarios x %d replicates x %.3f s/replicate (sustained, %d workers)",
    n_scen, R_FULL, sustained, WORKERS)
say("  PROJECTED WALL CLOCK: %.1f hours (%.1f days)", secs / 3600, secs / 86400)
say("")
bayes_share <- median(tm[, "bayes"], na.rm = TRUE) / tot
say("  Bayes is %.0f%% of serial cost. Two levers, both requiring evidence first:", 100 * bayes_share)
say("    - fitting ONE waning form instead of three cuts total time by ~%.0f%%,",
    100 * bayes_share * 2 / 3)
say("      but is only defensible if the pilot shows one form wins in >=95%% of")
say("      replicates. Pilot selection frequencies (1=constant, 2=affine in log t,")
say("      3=exponential decay):")
sel_all <- list()
for (nm in pilot_names) {
  dgn <- pilot[[nm]]$dgn
  b <- dgn[dgn$method == "bayes" & !is.na(dgn$bayes_form), ]
  if (!nrow(b)) { say("        %-16s no Bayesian fits recorded", nm); next }
  tb <- table(factor(b$bayes_form, levels = 1:3))
  fr <- 100 * tb / sum(tb)
  say("        %-16s form1 %4.0f%%  form2 %4.0f%%  form3 %4.0f%%   modal share %.0f%%",
      nm, fr[1], fr[2], fr[3], max(fr))
  say("        %-16s median |elpd margin| / SE = %.2f  (<1 means selection is",
      "", median(abs(b$bayes_elpd_margin_z), na.rm = TRUE))
  say("        %-16s  not distinguishable from a coin flip)", "")
  say("        %-16s replicates with any Pareto k > 0.7: %.0f%%", "",
      100 * mean(b$bayes_pareto_bad > 0, na.rm = TRUE))
  sel_all[[nm]] <- data.frame(scenario = nm, form = 1:3, pct = as.numeric(fr),
                              median_margin_z = median(abs(b$bayes_elpd_margin_z), na.rm = TRUE))
}
if (length(sel_all))
  write.csv(do.call(rbind, sel_all), "outputs/environment/bayes_form_selection_pilot.csv",
            row.names = FALSE)
say("")
say("  Tian bandwidth stability across the pilot:")
for (nm in pilot_names) {
  dgn <- pilot[[nm]]$dgn
  tb <- dgn[dgn$method == "tian" & !is.na(dgn$tian_bandwidth), ]
  if (!nrow(tb)) next
  say("        %-16s median h = %.1f d, IQR %.1f-%.1f, at search boundary in %.0f%% of replicates",
      nm, median(tb$tian_bandwidth), quantile(tb$tian_bandwidth, .25),
      quantile(tb$tian_bandwidth, .75), 100 * mean(tb$tian_h_at_boundary, na.rm = TRUE))
}
say("    - reusing adapted step size/metric across replicates, which must be")
say("      validated against full adaptation on a subsample before use.")
say("")
say("  Projection written to outputs/environment/runtime_projection.txt")
writeLines(c(
  sprintf("scenarios: %d", n_scen),
  sprintf("replicates: %d", R_FULL),
  sprintf("workers: %d", WORKERS),
  sprintf("sustained_sec_per_replicate: %.4f", sustained),
  sprintf("projected_hours: %.2f", secs / 3600),
  sprintf("bayes_share_of_serial_cost: %.3f", bayes_share)),
  "outputs/environment/runtime_projection.txt")
say("\n  DO NOT LAUNCH THE FULL RUN UNTIL THIS PROJECTION IS ACCEPTED (brief sec 5, Gate 4).")
