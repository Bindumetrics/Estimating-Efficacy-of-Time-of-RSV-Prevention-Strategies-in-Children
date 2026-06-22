# tests/test_summarise.R
# Deterministic shape/value checks for the summary tables on a synthetic
# per-day aggregate. Run from project root: Rscript tests/test_summarise.R

source(file.path("R", "evaluate.R"))
source(file.path("R", "summarise.R"))

# Minimal synthetic aggregated result: two methods, two days, one cell.
per_day <- expand.grid(scenario = "B", n = 2000L, loss_to_followup = 0,
                       method = c("M1", "M4"), day = c(90, 180),
                       stringsAsFactors = FALSE)
per_day$scale <- method_scale(per_day$method)
per_day$truth <- 0.7
per_day$mean_ve <- 0.72
per_day$sd_ve <- 0.05
per_day$mean_lower <- 0.6
per_day$mean_upper <- 0.85
per_day$signed_bias <- ifelse(per_day$method == "M1", 0.05, -0.03)
per_day$mae <- abs(per_day$signed_bias)
per_day$mse <- ifelse(per_day$method == "M1", 0.02, 0.01)
per_day$rmse <- sqrt(per_day$mse)
per_day$coverage <- 0.95
per_day$mean_width <- 0.25
per_day$n_replications_used <- 100L

out <- list(
  per_day = per_day,
  ise = data.frame(scenario = "B", n = 2000L, loss_to_followup = 0,
                   method = c("M1", "M4"), scale = c("hazard", "risk"),
                   mean_ise = c(2.0, 1.5), mean_iae = c(1.0, 0.8), n_ise = 100L),
  convergence = data.frame(scenario = "B", n = 2000L, loss_to_followup = 0,
                           method = c("M1", "M4"), n_replications = 100L,
                           convergence_rate = c(1.0, 0.98), mean_frac_days_available = c(1, 0.99)),
  m4_selection = data.frame(scenario = "B", n = 2000L, loss_to_followup = 0,
                            selected_form = c("exponential", "erlang3"), Freq = c(70L, 30L),
                            stringsAsFactors = FALSE)
)

b <- bias_by_day(out$per_day, days = c(90, 180))
stopifnot(all(c("d90", "d180") %in% names(b)))
stopifnot(nrow(b) == 2)                                    # two methods
stopifnot(b$d90[b$method == "M1"] == 0.05)

re <- relative_efficiency_table(out, reference = "M1", scale = "hazard")
stopifnot(re$rel_eff_vs_ref[re$method == "M1"] == 1)       # reference vs itself

cv <- convergence_table(out)
stopifnot(all(c("M1", "M4") %in% names(cv)))

sel <- m4_selection_table(out)
stopifnot(abs(sum(sel$proportion) - 1) < 1e-8)             # proportions sum to 1
stopifnot(abs(sel$proportion[sel$selected_form == "exponential"] - 0.7) < 1e-8)

cat("All summarise.R checks passed.\n")
