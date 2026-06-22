# run_full.R
# -----------------------------------------------------------------------------
# Full simulation pipeline driver: builds the FULL tiered design, runs the
# whole sweep (all five methods incl. M4/Stan), then writes every summary
# table and per-scenario figure. This is the real results run, not a smoke.
#
# Run from project root:  Rscript run_full.R
# (the run_full.bat launcher does this for you and finds Rscript on Windows.)
#
# Optional environment-variable overrides (all have sane defaults):
#   VE_PARALLEL   "1"/"true" to use multisession workers (default: TRUE)
#   VE_WORKERS    integer number of parallel workers (default: cores - 1)
#   VE_METHODS    comma list, e.g. "m1,m2,m3" to skip M4 (default: all five)
#   VE_FULL_R     override full-tier replications (default: config 300)
#   VE_HEAD_R     override headline-tier replications (default: config 1000)
#   VE_OUTDIR     sweep output dir (default: outputs/sweep)
#
# WARNING: the default full run is very large (~58k replications, M4/Stan
# dominates). Reduce VE_FULL_R / VE_HEAD_R, drop M4 via VE_METHODS, or raise
# VE_WORKERS for a tractable run. The sweep saves sweep_results.rds before the
# summary step, so summaries can be rebuilt without re-running the sweep.
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

# ---- Resolve options from environment --------------------------------------

env_flag <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || v == "") return(default)
  tolower(v) %in% c("1", "true", "yes", "y", "t")
}
env_int <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || v == "") return(default)
  as.integer(v)
}

parallel <- env_flag("VE_PARALLEL", TRUE)
workers  <- env_int("VE_WORKERS", max(1L, future::availableCores() - 1L))
out_dir  <- { v <- Sys.getenv("VE_OUTDIR", unset = ""); if (v == "") file.path("outputs", "sweep") else v }

methods_env <- Sys.getenv("VE_METHODS", unset = "")
methods <- if (methods_env == "") c("m1", "m1b", "m2", "m3", "m4") else
  trimws(strsplit(methods_env, ",")[[1]])

full_R <- env_int("VE_FULL_R", 300L)
head_R <- env_int("VE_HEAD_R", 1000L)

# ---- Build config, design, controls ----------------------------------------

config <- make_sim_config(m4_interval_width = 15L,
                          full_replications = full_R,
                          headline_replications = head_R)

design <- build_design(config)

# Under parallel multisession, run M4's chains sequentially inside each worker
# (parallel_chains = 1) to avoid oversubscription (workers x chains).
m4_ctrl <- list(chains = 4L,
                parallel_chains = if (parallel) 1L else 4L,
                iter_warmup = 1000L, iter_sampling = 1000L, adapt_delta = 0.95)
controls <- run_controls(methods = methods, m4 = m4_ctrl)

# ---- Report the plan before committing -------------------------------------

total_reps <- sum(design$n_replications)
cat("================ FULL SIMULATION RUN ================\n")
cat(sprintf("Cells:           %d  (full tier R=%d, headline tier R=%d)\n",
            nrow(design), full_R, head_R))
cat(sprintf("Total reps:      %d\n", total_reps))
cat(sprintf("Methods:         %s\n", paste(methods, collapse = ", ")))
cat(sprintf("Parallel:        %s%s\n", parallel,
            if (parallel) sprintf("  (workers = %d)", workers) else ""))
cat(sprintf("Sweep out dir:   %s\n", out_dir))
cat("=====================================================\n\n")

# ---- Run the sweep ----------------------------------------------------------

t0 <- Sys.time()
res <- run_simulation(design, config, controls = controls,
                      parallel = parallel, workers = workers,
                      out_dir = out_dir,
                      save_results = TRUE, save_raw_scored = FALSE, verbose = TRUE)
cat(sprintf("\nSweep elapsed: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- Summaries and figures --------------------------------------------------

summary_dir <- file.path("outputs", "summary")
figures_dir <- file.path("outputs", "figures")

tables <- write_summary_tables(res, out_dir = summary_dir)
figs   <- save_scenario_plots(res, out_dir = figures_dir)

cat("\nTables written to ", summary_dir, ":\n", sep = "")
cat("  ", paste(names(tables)[!vapply(tables, is.null, logical(1))], collapse = ", "), "\n")
cat("Figures written to ", figures_dir, ":\n", sep = "")
cat("  ", paste(basename(figs), collapse = ", "), "\n")

cat("\nFULL RUN COMPLETE.\n")
cat("  Sweep object : ", file.path(out_dir, "sweep_results.rds"), "\n", sep = "")
cat("  Tables       : ", summary_dir, "\n", sep = "")
cat("  Figures      : ", figures_dir, "\n", sep = "")
