# run_simulation.R
# -----------------------------------------------------------------------------
# Stage 4: tiered, parallel-safe orchestration of the simulation sweep.
#
# Design:
#   * One reproducible RNG stream per replication via L'Ecuyer-CMRG
#     (future.apply with future.seed = TRUE), seeded from config$meta$rng_seed.
#     Each replication also draws and STORES an integer `stan_seed` from its own
#     stream (used for M4's sampler and for record/audit).
#   * Tiered replications: the full factorial (all scenarios x n x censoring)
#     runs at `full_replications`; headline scenarios at the real-trial-sized
#     anchors run at `headline_replications`.
#   * Cell-by-cell execution bounds memory: a cell's replications are fitted,
#     scored, and aggregated, then the raw per-day scored data are discarded
#     (optionally saved). Only aggregates + small per-replication metadata are
#     retained, so a large sweep does not accumulate tens of millions of rows.
#   * Non-convergence rule: a method that fails (or errors) yields NA estimates
#     for that replication. NA estimates are excluded from point-estimate
#     aggregation but counted in the convergence-rate outcome (see evaluate.R).
# -----------------------------------------------------------------------------

# ---- Worker setup -----------------------------------------------------------

# Source all project code into a worker if it is not already loaded. Lets the
# same replication function run under sequential or multisession futures.
ensure_sources <- function(project_root) {
  if (exists("fit_m1", mode = "function") && exists("evaluate_methods", mode = "function") &&
      exists("simulate_trial", mode = "function")) {
    return(invisible())
  }
  source(file.path(project_root, "R", "config.R"))
  source(file.path(project_root, "R", "dgm.R"))
  source(file.path(project_root, "R", "evaluate.R"))
  for (m in c("m1", "m1b", "m2", "m3", "m4")) {
    source(file.path(project_root, "R", "methods", paste0(m, ".R")))
  }
  invisible()
}

# ---- Design grids -----------------------------------------------------------

# Full tiered design: every (scenario, n, censoring) cell, with headline cells
# (headline scenarios at the real-trial-sized anchors) upgraded to the headline
# replication count.
build_design <- function(config,
                         scenarios = scenario_ids(config),
                         sample_sizes = config$design$sample_sizes,
                         losses = config$design$loss_to_followup,
                         full_R = config$design$full_replications,
                         headline_R = config$design$headline_replications,
                         headline_scenarios = config$design$headline_scenarios,
                         headline_sizes = config$design$real_trial_size_anchors) {
  grid <- expand.grid(scenario = scenarios, n = as.integer(sample_sizes),
                      loss_to_followup = losses, stringsAsFactors = FALSE)
  is_headline <- grid$scenario %in% headline_scenarios & grid$n %in% headline_sizes
  grid$n_replications <- ifelse(is_headline, headline_R, full_R)
  grid$tier <- ifelse(is_headline, "headline", "full")
  grid[order(grid$scenario, grid$n, grid$loss_to_followup), ]
}

# Tiny design to validate the orchestration before any full-scale run.
smoke_design <- function(scenarios = c("B", "D"), sample_sizes = c(500L, 2000L),
                         losses = 0, n_replications = 3L) {
  grid <- expand.grid(scenario = scenarios, n = as.integer(sample_sizes),
                      loss_to_followup = losses, stringsAsFactors = FALSE)
  grid$n_replications <- as.integer(n_replications)
  grid$tier <- "smoke"
  grid
}

# ---- Run controls -----------------------------------------------------------

run_controls <- function(methods = c("m1", "m1b", "m2", "m3", "m4"),
                         m4 = list(chains = 4L, parallel_chains = 4L,
                                   iter_warmup = 1000L, iter_sampling = 1000L,
                                   adapt_delta = 0.95),
                         stan_file = file.path("stan", "m4_binomial.stan")) {
  list(methods = methods, m4 = m4, stan_file = stan_file)
}

# Lighter controls for smoke / CI runs.
smoke_controls <- function(methods = c("m1", "m1b", "m2", "m3", "m4")) {
  run_controls(methods = methods,
               m4 = list(chains = 2L, parallel_chains = 2L,
                         iter_warmup = 300L, iter_sampling = 300L, adapt_delta = 0.9))
}

# ---- One replication --------------------------------------------------------

safe_fit <- function(method_label, fn, ...) {
  res <- tryCatch(fn(...),
                  error = function(e) list(method = method_label, converged = FALSE,
                                           failure_reason = paste("unhandled error:", conditionMessage(e)),
                                           estimates = NULL))
  if (is.null(res)) {
    return(list(method = method_label, converged = FALSE,
                failure_reason = "fit returned NULL", estimates = NULL))
  }
  res
}

run_one_replication <- function(rep_id, scenario_id, calibrated, config, controls,
                                n, loss, model = NULL, project_root = getwd()) {
  grid <- seq.int(config$evaluation$grid_start, calibrated$spec$horizon)
  stan_seed <- sample.int(.Machine$integer.max, 1L)   # from this task's L'Ecuyer stream
  dat <- simulate_trial(n = n, calibrated = calibrated, loss_to_followup = loss, seed = NULL)

  M <- controls$methods
  fits <- list()
  if ("m1" %in% M)  fits$m1  <- safe_fit("M1", fit_m1, dat, grid = grid)
  if ("m1b" %in% M) fits$m1b <- safe_fit("M1b_supplementary", fit_m1b, dat, grid = grid)
  if ("m2" %in% M)  fits$m2  <- safe_fit("M2", fit_m2, dat, grid = grid,
                                         basis_dimension = config$m2$basis_dimension)
  if ("m3" %in% M)  fits$m3  <- safe_fit("M3", fit_m3, dat, grid = grid,
                                         bandwidth_days = config$m3$bandwidth_days,
                                         kernel = config$m3$kernel)
  if ("m4" %in% M) {
    if (is.null(model)) model <- m4_stan_model(controls$stan_file)
    fits$m4 <- safe_fit("M4", fit_m4, dat, grid = grid, horizon = calibrated$spec$horizon,
                        interval_width = config$m4$interval_width, model = model, seed = stan_seed,
                        chains = controls$m4$chains, parallel_chains = controls$m4$parallel_chains,
                        iter_warmup = controls$m4$iter_warmup,
                        iter_sampling = controls$m4$iter_sampling,
                        adapt_delta = controls$m4$adapt_delta, refresh = 0L)
  }

  ev <- evaluate_methods(fits, calibrated, scenario_id, config,
                         replication = rep_id, n = n, loss_to_followup = loss)

  meta <- data.frame(scenario = scenario_id, n = n, loss_to_followup = loss,
                     replication = rep_id, stan_seed = stan_seed,
                     m4_selected_form = if (!is.null(fits$m4)) fits$m4$selected_form %||% NA_character_ else NA_character_,
                     stringsAsFactors = FALSE)
  for (nm in names(fits)) meta[[paste0("conv_", fits[[nm]]$method)]] <- isTRUE(fits[[nm]]$converged)

  list(scored = ev$scored, ise = ev$ise, meta = meta)
}

# ---- M4 LOOIC-selection table ----------------------------------------------

# How often LOOIC selects each functional form, per cell (a result in itself:
# does LOOIC pick erlang3 under truth D, exponential under B, etc.?).
summarise_m4_selection <- function(meta) {
  if (!"m4_selected_form" %in% names(meta)) return(NULL)
  m <- meta[!is.na(meta$m4_selected_form), ]
  if (nrow(m) == 0) return(NULL)
  tab <- as.data.frame(table(scenario = m$scenario, n = m$n,
                             loss_to_followup = m$loss_to_followup,
                             selected_form = m$m4_selected_form),
                       stringsAsFactors = FALSE)
  tab <- tab[tab$Freq > 0, ]
  tab$n <- as.integer(tab$n)
  tab$loss_to_followup <- as.numeric(tab$loss_to_followup)
  tab[order(tab$scenario, tab$n, tab$loss_to_followup, tab$selected_form), ]
}

# ---- Top-level sweep --------------------------------------------------------

run_simulation <- function(design, config, controls = run_controls(),
                           parallel = FALSE,
                           workers = max(1L, future::availableCores() - 1L),
                           project_root = getwd(),
                           out_dir = file.path("outputs", "sweep"),
                           save_results = TRUE, save_raw_scored = FALSE,
                           verbose = TRUE) {
  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    stop("future and future.apply are required for the sweep.")
  }
  ensure_sources(project_root)
  calibrated_all <- calibrate_all_scenarios(config)

  # Compile the Stan model once up front (sequential plan reuses this object;
  # multisession workers re-load the cached executable themselves).
  model <- if ("m4" %in% controls$methods && !parallel) m4_stan_model(controls$stan_file) else NULL

  if (parallel) future::plan(future::multisession, workers = workers)
  else future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  set.seed(config$meta$rng_seed)   # makes the per-task L'Ecuyer streams reproducible

  per_day_list <- vector("list", nrow(design))
  conv_list <- vector("list", nrow(design))
  ise_list <- vector("list", nrow(design))
  meta_list <- vector("list", nrow(design))

  if (save_results) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(design))) {
    cell <- design[i, ]
    cal <- calibrated_all[[cell$scenario]]
    reps <- seq_len(cell$n_replications)
    if (verbose) {
      cat(sprintf("[%d/%d] scenario %s  n=%d  loss=%.2f  R=%d  tier=%s\n",
                  i, nrow(design), cell$scenario, cell$n, cell$loss_to_followup,
                  cell$n_replications, cell$tier))
    }

    rep_results <- future.apply::future_lapply(reps, function(r) {
      ensure_sources(project_root)
      run_one_replication(r, cell$scenario, cal, config, controls,
                          n = cell$n, loss = cell$loss_to_followup,
                          model = model, project_root = project_root)
    }, future.seed = TRUE)

    scored <- do.call(rbind, lapply(rep_results, `[[`, "scored"))
    ise <- do.call(rbind, lapply(rep_results, `[[`, "ise"))
    meta <- do.call(rbind, lapply(rep_results, `[[`, "meta"))

    per_day_list[[i]] <- aggregate_per_day(scored)
    conv_list[[i]] <- aggregate_convergence(scored)
    ise_list[[i]] <- aggregate_ise(ise)
    meta_list[[i]] <- meta

    if (save_raw_scored) {
      saveRDS(scored, file.path(out_dir, sprintf("scored_%s_n%d_l%02d.rds",
                                                 cell$scenario, cell$n,
                                                 round(cell$loss_to_followup * 100))))
    }
  }

  out <- list(
    design = design,
    per_day = do.call(rbind, per_day_list),
    convergence = do.call(rbind, conv_list),
    ise = do.call(rbind, ise_list),
    meta = do.call(rbind, meta_list),
    config = config,
    controls = controls
  )
  out$m4_selection <- summarise_m4_selection(out$meta)

  if (save_results) {
    saveRDS(out, file.path(out_dir, "sweep_results.rds"))
    if (verbose) cat("Saved sweep results to", file.path(out_dir, "sweep_results.rds"), "\n")
  }
  out
}
