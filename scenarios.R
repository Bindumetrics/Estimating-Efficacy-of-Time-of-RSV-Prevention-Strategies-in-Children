# ---------------------------------------------------------------------------
# The scenario grid: one factor varied at a time around a reference.
#
# Reference (brief sec 3.3, revised): MATISSE FINAL analysis -- arms 3585/3563,
# 180-day horizon, ~200 events, 4% LTFU plus staggered-entry administrative
# censoring, gradual waning, 15-day evaluation grid, MATISSE-derived seasonal
# baseline hazard.
#
# Every entry carries a `tag` of "verified" or "assumed". Anything tagged
# assumed must appear in the limitations section of REPORT.md.
# ---------------------------------------------------------------------------

build_scenarios <- function() {
  S <- list()
  add <- function(...) { s <- make_scenario(...); S[[s$name]] <<- s; invisible(NULL) }

  # --- reference -----------------------------------------------------------
  add("ref", notes = "verified: arms, horizon, events, LTFU, baseline hazard")

  # --- factor: VE shape ----------------------------------------------------
  for (sh in c("constant", "rapid", "piecewise")) {
    add(paste0("shape_", sh), ve_shape = sh,
        notes = "assumed: VE shape chosen for design coverage, not from data")
  }
  # frailty-induced: VE_id constant, target is the marginalised VE_h
  for (th in c(0.5, 1.0)) {
    add(sprintf("frailty_%.1f", th), ve_shape = "constant", theta = th,
        notes = "assumed: gamma frailty variance; produces apparent waning with no individual waning")
  }

  # --- factor: event count -------------------------------------------------
  for (ev in c(30, 70, 93, 304)) {
    add(paste0("events_", ev), target_events = ev,
        notes = if (ev %in% c(93, 304)) "verified: MATISSE severe @180d / any-severity @360d"
                else "assumed: sparse/low event counts for stress testing")
  }

  # --- factor: sample size -------------------------------------------------
  add("size_nirsevimab", n_c = 786, n_v = 1564, horizon = 150,
      target_events = 70, ltfu_prob = 0.02, accrual_days = 90,
      notes = "verified: Simoes 2023 pooled -- arms, horizon, events, 2% dropout")
  add("size_doubled", n_c = 7126, n_v = 7170,
      notes = "assumed: doubled arms, synthetic")

  # --- factor: loss to follow-up, crossed with administrative censoring -----
  for (l in c(0, 0.10, 0.20)) {
    add(sprintf("ltfu_%02d", round(100 * l)), ltfu_prob = l,
        notes = "assumed: LTFU sweep; scale held at reference (see below)")
  }
  add("admin_off", accrual_days = 0,
      notes = "assumed: administrative censoring disabled -- the configuration that degenerates Tian's kernel weights")
  add("admin_off_noltfu", accrual_days = 0, ltfu_prob = 0,
      notes = "assumed: DEGENERATE case, all non-cases exit at the horizon")

  # --- factor: horizon -----------------------------------------------------
  add("horizon_150", horizon = 150, target_events = 164,
      notes = "verified: MATISSE any-severity @150d")
  add("horizon_360", horizon = 360, target_events = 304, bayes_interval = 45,
      pl_interval = 60,
      notes = "verified: MATISSE any-severity @360d -- NOT synthetic")

  # --- factor: baseline hazard shape ---------------------------------------
  add("hazard_flat", h0 = h0_flat, notes = "assumed: flat baseline")
  add("hazard_increasing", h0 = h0_increasing, notes = "assumed: monotone increasing baseline")

  # --- applied anchors -----------------------------------------------------
  fitobj <- tryCatch(readRDS("outputs/environment/applied_anchor_fit.rds"), error = function(e) NULL)
  if (!is.null(fitobj)) {
    fitfn <- ve_logquad(fitobj$par[1], fitobj$par[2], fitobj$par[3])
    add("applied_matisse_180", ve_shape = "applied", ve_fn = fitfn,
        target_events = 200, horizon = 180,
        notes = "verified: VE_id calibrated to the 8-point published cumulative VE trajectory")
    add("applied_matisse_360", ve_shape = "applied", ve_fn = fitfn,
        target_events = 304, horizon = 360, bayes_interval = 45, pl_interval = 60,
        notes = "verified: as above, 360-day horizon")
  }
  S
}

# The LTFU sweep must hold the baseline scale FIXED at the reference level so
# that higher dropout genuinely yields fewer events, which is the question the
# sweep claims to answer (brief sec 3.1). Every other sweep re-solves the scale.
fix_ltfu_sweep_scale <- function(S) {
  ref_scale <- S[["ref"]]$scale
  for (nm in grep("^ltfu_|^admin_off", names(S), value = TRUE)) {
    S[[nm]]$scale <- ref_scale
    S[[nm]]$h0_scaled <- list(cuts = S[[nm]]$h0$cuts, vals = S[[nm]]$h0$vals * ref_scale)
    S[[nm]]$ve_h_true <- ve_h_true(S[[nm]]$grid, S[[nm]]$h0_scaled,
                                   S[[nm]]$ve_fn, S[[nm]]$theta)
    S[[nm]]$notes <- paste(S[[nm]]$notes, "| scale HELD FIXED at reference")
  }
  S
}

scenario_register <- function(S) {
  do.call(rbind, lapply(S, function(s) data.frame(
    scenario = s$name, n_c = s$n_c, n_v = s$n_v, horizon = s$horizon,
    target_events = s$target_events, ve_shape = s$ve_shape, theta = s$theta,
    ltfu_prob = s$ltfu_prob, accrual_days = s$accrual_days,
    baseline_scale = s$scale, grid_points = length(s$grid),
    pl_interval = s$pl_interval, bayes_interval = s$bayes_interval,
    tag = ifelse(grepl("^verified", s$notes), "verified", "assumed"),
    notes = s$notes, stringsAsFactors = FALSE)))
}
