# tests/test_evaluate.R
# Deterministic checks of the evaluation metric math (no Stan / no fitting).
# Run from project root: Rscript tests/test_evaluate.R

source(file.path("R", "config.R"))
source(file.path("R", "dgm.R"))
source(file.path("R", "evaluate.R"))

approx_eq <- function(a, b, tol = 1e-8) all(abs(a - b) < tol)

# --- method_scale ------------------------------------------------------------
stopifnot(method_scale("M1") == "hazard")
stopifnot(method_scale("M2") == "hazard")
stopifnot(method_scale("M3") == "hazard")
stopifnot(method_scale("M1b_supplementary") == "hazard")
stopifnot(method_scale("M4") == "risk")

# --- day_to_interval ---------------------------------------------------------
ig <- m4_interval_grid(60, 15L)           # intervals (0,15],(15,30],(30,45],(45,60]
stopifnot(day_to_interval(c(1, 15, 16, 30, 45, 60), ig) == c(1, 1, 2, 2, 3, 4))

# --- trapz -------------------------------------------------------------------
stopifnot(approx_eq(trapz(1:5, rep(2, 5)), 2 * 4))         # constant 2 over width 4
stopifnot(approx_eq(trapz(c(0, 2), c(0, 2)), 2))            # triangle area
stopifnot(is.na(trapz(1, 1)))                              # < 2 points

# --- harmonise + attach_truth + score ---------------------------------------
grid <- 1:4
targets <- data.frame(day = grid,
                      truth_hazard = c(0.90, 0.80, 0.70, 0.60),
                      truth_risk   = c(0.95, 0.85, 0.75, 0.65))
fit_m1 <- list(method = "M1", converged = TRUE,
               estimates = data.frame(day = grid, ve = c(0.92, 0.78, 0.72, 0.55),
                                      lower = c(0.82, 0.68, 0.62, 0.45),
                                      upper = c(1.02, 0.88, 0.82, 0.65)))
fit_m4 <- list(method = "M4", converged = TRUE,
               estimates = data.frame(day = grid, ve = c(0.90, 0.80, 0.78, 0.70),
                                      lower = c(0.85, 0.75, 0.73, 0.65),
                                      upper = c(0.95, 0.85, 0.83, 0.75)))
fit_m2 <- list(method = "M2", converged = FALSE, estimates = NULL)   # total failure

h <- harmonise_estimates(list(fit_m1, fit_m4, fit_m2), grid)
stopifnot(nrow(h) == 12)
stopifnot(all(is.na(h$ve[h$method == "M2"])))

scored <- score_estimates(attach_truth(h, targets))
m1 <- scored[scored$method == "M1", ]
stopifnot(approx_eq(m1$error, c(0.02, -0.02, 0.02, -0.05)))   # vs hazard truth
stopifnot(all(m1$covered == 1L))
stopifnot(approx_eq(m1$width, rep(0.20, 4)))

m4 <- scored[scored$method == "M4", ]
stopifnot(approx_eq(m4$error, c(-0.05, -0.05, 0.03, 0.05)))   # vs risk truth (routing works)
stopifnot(all(m4$covered == 1L))

m2 <- scored[scored$method == "M2", ]
stopifnot(all(!m2$available), all(is.na(m2$covered)))

# --- replication_ise ---------------------------------------------------------
ise <- replication_ise(scored)
m1_ise <- ise$ise[ise$method == "M1"]
expected_m1_ise <- trapz(grid, m1$sq_error)
stopifnot(approx_eq(m1_ise, expected_m1_ise))
stopifnot(is.na(ise$ise[ise$method == "M2"]))                  # no available days

# --- aggregation across replications ----------------------------------------
mk_scored <- function(rep_id, errs) {
  data.frame(scenario = "B", n = 500L, loss_to_followup = 0,
             replication = rep_id, method = "M1", scale = "hazard", day = grid,
             ve = NA, lower = NA, upper = NA, truth = 0.8, truth_common = 0.8,
             available = TRUE, error = errs, abs_error = abs(errs),
             sq_error = errs^2, width = 0.2, covered = 1L,
             error_c = errs, abs_error_c = abs(errs), sq_error_c = errs^2, covered_c = 1L)
}
stack2 <- rbind(mk_scored(1L, c(0.1, 0.1, 0.1, 0.1)),
                mk_scored(2L, c(-0.1, -0.1, -0.1, -0.1)))
agg <- aggregate_per_day(stack2)
stopifnot(approx_eq(agg$signed_bias, rep(0, 4)))               # +0.1 and -0.1 cancel
stopifnot(approx_eq(agg$mse, rep(0.01, 4)))                    # both 0.1^2
stopifnot(approx_eq(agg$coverage, rep(1, 4)))
stopifnot(all(agg$n_replications_used == 2))

# convergence: add a failed replication (no available days)
failed_rep <- mk_scored(3L, rep(0, 4)); failed_rep$available <- FALSE
conv <- aggregate_convergence(rbind(stack2, failed_rep))
stopifnot(approx_eq(conv$convergence_rate, 2 / 3))             # 2 of 3 reps usable

# relative efficiency: M_a MSE twice M_b -> ratio 2
re_in <- rbind(
  data.frame(scenario = "B", n = 500L, loss_to_followup = 0, method = "M_a",
             scale = "hazard", day = grid, mse = 0.02),
  data.frame(scenario = "B", n = 500L, loss_to_followup = 0, method = "M_b",
             scale = "hazard", day = grid, mse = 0.01)
)
re <- relative_efficiency(re_in, "B", 500L, 0, scale = "hazard")
stopifnot(approx_eq(re$ratio["M_a", "M_b"], 2))
stopifnot(approx_eq(re$ratio["M_b", "M_a"], 0.5))

# --- common-scale scoring ----------------------------------------------------
# A risk-scale (M4-like) method scored on the common hazard reference must use
# truth_hazard, not its own risk truth.
hM4 <- harmonise_estimates(list(fit_m4), grid)
sc4 <- score_estimates(attach_truth(hM4, targets, common_scale = "hazard"))
stopifnot(approx_eq(sc4$truth_common, targets$truth_hazard))
stopifnot(approx_eq(sc4$error_c, fit_m4$estimates$ve - targets$truth_hazard))
stopifnot(approx_eq(sc4$error, fit_m4$estimates$ve - targets$truth_risk))   # own scale unchanged

# --- M4 per-form expansion ---------------------------------------------------
m4_fit <- list(
  method = "M4", converged = TRUE,
  estimates = data.frame(day = grid, ve = 0.7, lower = 0.6, upper = 0.8),
  selected_form = "exponential",
  per_form = list(
    exponential = list(converged = TRUE,
                       estimates = data.frame(method = "M4", form = "exponential",
                                              day = grid, ve = 0.71, lower = 0.6, upper = 0.8)),
    erlang3 = list(converged = TRUE,
                   estimates = data.frame(method = "M4", form = "erlang3",
                                          day = grid, ve = 0.65, lower = 0.5, upper = 0.78)),
    power_law = list(converged = FALSE, estimates = NULL)
  )
)
expanded <- expand_fits_for_scoring(list(m4 = m4_fit), score_m4_per_form = TRUE)
stopifnot(all(c("M4", "M4_exponential", "M4_erlang3", "M4_power_law") %in% names(expanded)))
stopifnot(method_scale("M4_exponential") == "risk")
stopifnot(!isTRUE(expanded$M4_power_law$converged))               # failed form preserved
# With the flag off, the M4 fit passes through unexpanded.
stopifnot(identical(names(expand_fits_for_scoring(list(m4 = m4_fit), FALSE)), "m4"))

# --- truth_targets integration (real DGM) -----------------------------------
config <- make_sim_config()
calibrated <- calibrate_all_scenarios(config)
tt <- truth_targets(calibrated[["B"]], "B", config)
stopifnot(nrow(tt) == calibrated[["B"]]$spec$horizon)
stopifnot(all(is.finite(tt$truth_hazard)), all(is.finite(tt$truth_risk)))
# Under rare events the risk-scale (interval-conditional) VE sits above the
# hazard-scale VE; both should be ordered sensibly and within (0, 1).
stopifnot(all(tt$truth_hazard > 0 & tt$truth_hazard < 1))
stopifnot(all(tt$truth_risk > 0 & tt$truth_risk < 1))

# Scale gap should be small (rare events) and finite for every scenario.
gap <- scale_gap_table(config)
stopifnot(all(is.finite(gap$max_abs_gap)))
stopifnot(all(gap$max_abs_gap < 0.25))                            # scales stay close

cat("All evaluate.R checks passed.\n")
