# validate_summary.R
# Stage 5 validation: build all summary tables and per-scenario figures from a
# saved sweep result (the Stage 4 smoke run by default).
#
# Run from project root:  Rscript validate_summary.R

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))
source(file.path("R", "evaluate.R"))
source(file.path("R", "summarise.R"))
source(file.path("R", "plots.R"))

sweep_file <- file.path("outputs", "sweep_smoke", "sweep_results.rds")
if (!file.exists(sweep_file)) stop("Run validate_run.R first to produce ", sweep_file)
out <- readRDS(sweep_file)

tables <- write_summary_tables(out, out_dir = file.path("outputs", "summary_smoke"))
figs <- save_scenario_plots(out, out_dir = file.path("outputs", "figures_smoke"))

cat("Tables written to outputs/summary_smoke:\n")
cat(" ", paste(names(tables)[!vapply(tables, is.null, logical(1))], collapse = ", "), "\n\n")

cat("Signed bias by day (direction + magnitude):\n")
print(tables$bias, row.names = FALSE)

cat("\n95% coverage by day:\n")
print(tables$coverage, row.names = FALSE)

cat("\nScale gap between VE_h and VE_r (justifies the common-scale comparison):\n")
print(tables$scale_gap, row.names = FALSE)

cat("\nFOUR-WAY RANKING (common hazard scale, all methods; rank 1 = most efficient):\n")
rec <- tables$relative_efficiency_common
print(rec[, c("scenario", "n", "method", "integrated_mse", "rel_eff_vs_ref", "rank")], row.names = FALSE)

cat("\nRelative efficiency vs M1 (own hazard scale, integrated MSE ratio; >1 = worse):\n")
re <- tables$relative_efficiency
print(re[, c("scenario", "n", "method", "integrated_mse", "rel_eff_vs_ref", "rank")], row.names = FALSE)

cat("\nM4 per-form signed bias by day (misspecification tier; each form vs every truth):\n")
print(tables$m4_per_form_bias, row.names = FALSE)

cat("\nM4 LOOIC form-selection proportions:\n")
print(tables$m4_selection[, c("scenario", "n", "selected_form", "Freq", "proportion")], row.names = FALSE)

cat("\nFigures written:\n")
cat(" ", paste(basename(figs), collapse = ", "), "\n")
cat("\nStage 5 validation complete. Tables in outputs/summary_smoke, figures in outputs/figures_smoke.\n")
