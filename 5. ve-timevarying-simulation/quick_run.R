# quick_run.R
# -----------------------------------------------------------------------------
# Deadline-friendly reduced run: same sweep + summary + figure code as
# run_full.R, but with a capped sample-size ladder and low replications so the
# whole thing finishes in minutes instead of hours. Drops M4/Stan by default.
#
#   Rscript quick_run.R
#
# This does NOT touch run_full.R or config.R. It writes the sweep to
# outputs/sweep_quick and the summaries/figures to the standard
# outputs/summary and outputs/figures (so a report reading those picks them up).
# -----------------------------------------------------------------------------

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))
source(file.path("R", "evaluate.R"))
source(file.path("R", "methods", "m1.R"))
source(file.path("R", "methods", "m1b.R"))
source(file.path("R", "methods", "m2.R"))
source(file.path("R", "methods", "m3.R"))
source(file.path("R", "methods", "m4.R"))
source(file.path("R", "run_simulation.R"))
source(file.path("R", "summarise.R"))
source(file.path("R", "plots.R"))

# ---- Reduced knobs ----------------------------------------------------------
quick_sample_sizes <- c(500L, 2000L, 5000L)   # drops the costly 7000 & 10000
quick_R            <- 10L                       # low replications for speed
quick_methods      <- c("m1", "m1b", "m2", "m3")  # no M4/Stan
parallel <- TRUE
workers  <- max(1L, future::availableCores() - 1L)
out_dir  <- file.path("outputs", "sweep_quick")

config <- make_sim_config(full_replications = quick_R, headline_replications = quick_R)

# Cap the sample-size ladder; with these sizes nothing hits the headline
# anchors (2350/7000), so every cell runs at quick_R in the "full" tier.
design <- build_design(config, sample_sizes = quick_sample_sizes,
                       full_R = quick_R, headline_R = quick_R)

controls <- run_controls(methods = quick_methods)

total_reps <- sum(design$n_replications)
cat("================ QUICK SIMULATION RUN ================\n")
cat(sprintf("Cells:           %d  (R=%d per cell)\n", nrow(design), quick_R))
cat(sprintf("Sample sizes:    %s\n", paste(quick_sample_sizes, collapse = ", ")))
cat(sprintf("Total reps:      %d\n", total_reps))
cat(sprintf("Methods:         %s\n", paste(quick_methods, collapse = ", ")))
cat(sprintf("Parallel:        %s  (workers = %d)\n", parallel, workers))
cat(sprintf("Sweep out dir:   %s\n", out_dir))
cat("=====================================================\n\n")

t0 <- Sys.time()
res <- run_simulation(design, config, controls = controls,
                      parallel = parallel, workers = workers,
                      out_dir = out_dir,
                      save_results = TRUE, save_raw_scored = FALSE, verbose = TRUE)
cat(sprintf("\nSweep elapsed: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

summary_dir <- file.path("outputs", "summary")
figures_dir <- file.path("outputs", "figures")
tables <- write_summary_tables(res, out_dir = summary_dir)
figs   <- save_scenario_plots(res, out_dir = figures_dir)

cat("\nTables written to ", summary_dir, "\n", sep = "")
cat("Figures written to ", figures_dir, "\n", sep = "")
cat("\nQUICK RUN COMPLETE.\n")
