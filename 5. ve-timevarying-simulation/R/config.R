# config.R
# -----------------------------------------------------------------------------
# Stage 1 configuration for the thesis simulation study.
#
# Statistical aim: estimate how well each method recovers a known time-varying
# vaccine effectiveness curve. This is not a power simulation.
# -----------------------------------------------------------------------------

make_sim_config <- function(
  m4_interval_width = 15L,
  full_replications = 300L,
  headline_replications = 1000L,
  rng_seed = 20260613L
) {
  stopifnot(m4_interval_width > 0)

  list(
    meta = list(
      project = "time-varying VE simulation study",
      stage = "stage 1: config + data-generating mechanism",
      created = as.character(Sys.Date()),
      rng_seed = rng_seed
    ),

    evaluation = list(
      grid_start = 1L,
      summary_days = c(30L, 60L, 90L, 120L, 150L, 180L),
      compare_m4_on_risk_scale = TRUE,
      note = paste(
        "M1-M3 target hazard-scale VE_h(t). M4 targets risk-scale VE_r(t).",
        "The DGM stores both scales so evaluation can avoid confounding",
        "method bias with a pure scale difference."
      )
    ),

    design = list(
      sample_sizes = c(500L, 2000L, 2350L, 5000L, 7000L, 10000L),
      real_trial_size_anchors = c(2350L, 7000L),
      loss_to_followup = c(0, 0.05, 0.10, 0.20),
      full_replications = full_replications,
      headline_replications = headline_replications,
      headline_scenarios = c("A", "B", "C", "D")
    ),

    m2 = list(
      smoother = "mgcv penalised spline on cox.zph scaled Schoenfeld estimates",
      basis_dimension = 8L,
      method = "REML",
      ci = "pointwise 95% intervals from smoother standard errors on the log-HR scale"
    ),

    m3 = list(
      estimator = "Tian / Zucker-Wei kernel-weighted local partial likelihood",
      kernel = "gaussian",
      bandwidth_days = 30,
      ci = "pointwise 95% Wald intervals from the Tian/Zucker-Wei sandwich variance A^{-1} B A^{-1} (squared-kernel meat)",
      note = "No turnkey CRAN package is used; this is implemented directly. Bandwidth sensitivity belongs outside the primary run."
    ),

    m4 = list(
      interval_width = as.integer(m4_interval_width),
      candidate_forms = c("exponential", "erlang3", "power_law"),
      primary_prior = list(
        ve0 = "Beta(2, 2)",
        kappa = "Half-Normal, weakly informative",
        erlang_rate = "Gamma(2, 1)",
        note = "Matched/informative Beta(8, 2) is reserved for sensitivity analyses."
      ),
      likelihood = "binomial interval-count likelihood with pi_vk = pi_pk * (1 - VE_r,k)"
    ),

    baseline = list(
      hazard_family = "Weibull",
      shape_k = 1.3,
      lambda_rule = "lambda is calibrated from target placebo cumulative incidence: -log(1 - p) / horizon^k"
    ),

    trial_anchors = list(
      nirsevimab = list(
        source = "Simões et al. (2023), pooled nirsevimab; thesis Ch. 3 / Section 5.2 anchor",
        cumulative_ve = 0.795,
        horizon = 150L,
        placebo_cumulative_incidence = 0.060,
        randomisation_ratio = c(vaccine = 2L, placebo = 1L),
        note = paste(
          "Scenario A is a specificity/null-control scenario for false waning.",
          "It is not a claim that nirsevimab VE is constant."
        )
      ),
      matisse = list(
        source = "Kampmann et al. (2023), MATISSE severe RSV-LRTI cumulative VE anchors",
        cumulative_ve = c(`90` = 0.818, `180` = 0.694),
        horizon = 180L,
        placebo_cumulative_incidence = 0.030,
        randomisation_ratio = c(vaccine = 1L, placebo = 1L),
        note = paste(
          "Only the 90-day and 180-day severe RSV-LRTI anchors are used by default.",
          "Intermediate anchors should be added only if directly verified from the source."
        )
      )
    ),

    scenarios = list(
      A = list(
        label = "Constant VE specificity control",
        family = "constant",
        horizon = 150L,
        target_placebo_cumulative_incidence = 0.060,
        randomisation_ratio = c(vaccine = 2L, placebo = 1L),
        fixed = list(ve = 0.795),
        anchors = c(`150` = 0.795),
        role = "Specificity control: tests whether methods falsely infer waning."
      ),
      B = list(
        label = "Exponential waning, MATISSE magnitude anchor",
        family = "exponential",
        horizon = 180L,
        target_placebo_cumulative_incidence = 0.030,
        randomisation_ratio = c(vaccine = 1L, placebo = 1L),
        anchors = c(`90` = 0.818, `180` = 0.694),
        initial = list(ve0 = 0.90, rho = 0.0025),
        role = "Well-specified truth for M4 exponential candidate."
      ),
      C = list(
        label = "Weibull/power-law-like waning, MATISSE magnitude anchor",
        family = "weibull",
        horizon = 180L,
        target_placebo_cumulative_incidence = 0.030,
        randomisation_ratio = c(vaccine = 1L, placebo = 1L),
        anchors = c(`90` = 0.818, `180` = 0.694),
        fixed = list(gamma = 1.5),
        initial = list(ve0 = 0.90, rho = 0.0040),
        role = "Near-but-not-exact truth for M4 power-law candidate."
      ),
      D = list(
        label = "Erlang-3 / biphasic maternal-antibody-like waning",
        family = "erlang3",
        horizon = 180L,
        target_placebo_cumulative_incidence = 0.030,
        randomisation_ratio = c(vaccine = 1L, placebo = 1L),
        anchors = c(`90` = 0.818, `180` = 0.694),
        initial = list(ve0 = 0.90, rate = 0.020),
        role = "Well-specified truth for M4 Erlang-3 candidate."
      ),
      E = list(
        label = "Delayed waning stress test",
        family = "delayed_exponential",
        horizon = 180L,
        target_placebo_cumulative_incidence = 0.030,
        randomisation_ratio = c(vaccine = 1L, placebo = 1L),
        anchors = c(`90` = 0.818, `180` = 0.694),
        fixed = list(delay = 120),
        initial = list(ve0 = 0.82, rho = 0.010),
        role = "Methodological stress test; no direct trial-shape anchor."
      )
    )
  )
}

scenario_ids <- function(config) {
  names(config$scenarios)
}

m4_interval_grid <- function(horizon, interval_width = 15L) {
  breaks <- unique(c(0L, seq.int(interval_width, horizon, by = interval_width), horizon))
  breaks <- sort(breaks)
  data.frame(
    interval = seq_len(length(breaks) - 1L),
    start = head(breaks, -1L),
    end = tail(breaks, -1L)
  )
}


