# tests/test_run.R
# Fast orchestration checks: design tiering, sweep plumbing, non-convergence
# handling, and reproducibility. M4 is excluded so the test needs no Stan.
# Run from project root: Rscript tests/test_run.R

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))
source(file.path("R", "evaluate.R"))
source(file.path("R", "methods", "m1.R"))
source(file.path("R", "methods", "m1b.R"))
source(file.path("R", "methods", "m2.R"))
source(file.path("R", "methods", "m3.R"))
source(file.path("R", "methods", "m4.R"))
source(file.path("R", "run_simulation.R"))

config <- make_sim_config()

# --- design tiering ----------------------------------------------------------
design <- build_design(config)
stopifnot(nrow(design) == length(scenario_ids(config)) *
            length(config$design$sample_sizes) * length(config$design$loss_to_followup))
headline <- design[design$scenario == "B" & design$n == 7000L, ]
full <- design[design$scenario == "E" & design$n == 500L, ]
stopifnot(all(headline$n_replications == config$design$headline_replications))
stopifnot(all(full$n_replications == config$design$full_replications))
stopifnot(all(design$tier[design$scenario == "E"] == "full"))   # E never headline
cat("Design tiering checks passed.\n")

# --- small sweep (no Stan) ---------------------------------------------------
d <- smoke_design(scenarios = "B", sample_sizes = 2000L, losses = 0, n_replications = 2L)
ctrl <- run_controls(methods = c("m1", "m2", "m3"))

res1 <- run_simulation(d, config, controls = ctrl, parallel = FALSE,
                       save_results = FALSE, verbose = FALSE)

stopifnot(all(c("per_day", "convergence", "ise", "meta", "design") %in% names(res1)))
stopifnot(nrow(res1$per_day) > 0)
stopifnot(all(c("signed_bias", "mse", "coverage", "mean_width", "mean_ve") %in% names(res1$per_day)))
stopifnot(all(res1$convergence$convergence_rate >= 0 & res1$convergence$convergence_rate <= 1))
stopifnot(all(c("M1", "M2", "M3") %in% res1$per_day$method))
stopifnot(!any(res1$per_day$method == "M4"))                     # M4 excluded
stopifnot("stan_seed" %in% names(res1$meta))                    # per-replication seed stored
stopifnot(nrow(res1$meta) == 2L)                               # 2 replications
cat("Sweep plumbing checks passed.\n")

# --- reproducibility ---------------------------------------------------------
res2 <- run_simulation(d, config, controls = ctrl, parallel = FALSE,
                       save_results = FALSE, verbose = FALSE)
m1 <- res1$per_day[res1$per_day$method == "M1", ]
m1b <- res2$per_day[res2$per_day$method == "M1", ]
m1 <- m1[order(m1$day), ]; m1b <- m1b[order(m1b$day), ]
stopifnot(isTRUE(all.equal(m1$signed_bias, m1b$signed_bias)))
stopifnot(isTRUE(all.equal(res1$meta$stan_seed, res2$meta$stan_seed)))
cat("Reproducibility checks passed.\n")

cat("All run_simulation.R checks passed.\n")
