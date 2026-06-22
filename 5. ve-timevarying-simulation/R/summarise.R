# summarise.R
# -----------------------------------------------------------------------------
# Stage 5 (tables): presentation tables from a run_simulation() result.
#
# Bias tables show DIRECTION and MAGNITUDE per method x scenario x n (not just
# RMSE magnitudes). Separate tables give RMSE, coverage, and mean interval width
# by the standard reporting days; ISE; pairwise relative efficiency (within a
# truth scale, since M4's risk-scale estimand is not directly comparable to the
# hazard-scale methods); convergence rate; and the M4 LOOIC functional-form
# selection table. Everything is written as CSV; headline tables are also
# rendered as gt HTML.
# -----------------------------------------------------------------------------

# Wide metric-by-day table: one row per (scenario, n, loss, method, scale),
# one column per reporting day.
metric_by_day <- function(per_day, metric, days = c(30, 60, 90, 120, 150, 180)) {
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr is required.")
  days <- days[days %in% per_day$day]
  cols <- c("scenario", "n", "loss_to_followup", "method", "scale", "day", metric)
  sub <- per_day[per_day$day %in% days, cols]
  names(sub)[names(sub) == metric] <- "value"
  out <- tidyr::pivot_wider(sub, names_from = "day", values_from = "value",
                            names_prefix = "d")
  out <- out[order(out$scenario, out$n, out$loss_to_followup, out$method), ]
  as.data.frame(out)
}

bias_by_day     <- function(per_day, days = c(30, 60, 90, 120, 150, 180)) metric_by_day(per_day, "signed_bias", days)
rmse_by_day     <- function(per_day, days = c(30, 60, 90, 120, 150, 180)) metric_by_day(per_day, "rmse", days)
coverage_by_day <- function(per_day, days = c(30, 60, 90, 120, 150, 180)) metric_by_day(per_day, "coverage", days)
width_by_day    <- function(per_day, days = c(30, 60, 90, 120, 150, 180)) metric_by_day(per_day, "mean_width", days)

# ISE table (mean integrated squared error per cell x method). M4 is on the risk
# scale; the column is flagged so it is not read as directly comparable.
ise_table <- function(out) {
  t <- out$ise[, c("scenario", "n", "loss_to_followup", "method", "scale", "mean_ise", "mean_iae", "n_ise")]
  t[order(t$scenario, t$n, t$loss_to_followup, t$method), ]
}

# Pairwise relative efficiency: horizon-integrated MSE per method, the ratio to a
# reference method, and the within-cell rank (1 = lowest MSE = most efficient).
#   * own-scale  (scale = "hazard", metric = "mse"): bias-pure, hazard methods.
#   * common-scale (scale = NULL, metric = "mse_common"): all four methods ranked
#     head-to-head. This is the RQ4 ranking; pair with scale_gap_table().
relative_efficiency_table <- function(out, reference = "M1", scale = "hazard",
                                      metric = "mse", include_m4_forms = FALSE) {
  pd <- out$per_day
  if (!is.null(scale)) pd <- pd[pd$scale == scale, ]
  if (!include_m4_forms) pd <- pd[!grepl("^M4_", pd$method), ]
  if (nrow(pd) == 0) return(NULL)
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required.")
  pd$.metric <- pd[[metric]]
  imse <- dplyr::summarise(
    dplyr::group_by(pd, .data$scenario, .data$n, .data$loss_to_followup, .data$method),
    integrated_mse = mean(.data$.metric), .groups = "drop"
  )
  imse <- as.data.frame(imse)
  ref <- imse[imse$method == reference, c("scenario", "n", "loss_to_followup", "integrated_mse")]
  names(ref)[names(ref) == "integrated_mse"] <- "ref_mse"
  m <- merge(imse, ref, by = c("scenario", "n", "loss_to_followup"), all.x = TRUE)
  m$rel_eff_vs_ref <- m$integrated_mse / m$ref_mse        # > 1 means worse than reference
  m$reference <- reference
  m <- do.call(rbind, by(m, list(m$scenario, m$n, m$loss_to_followup), function(d) {
    d$rank <- rank(d$integrated_mse, ties.method = "min")
    d
  }))
  m[order(m$scenario, m$n, m$loss_to_followup, m$rank), ]
}

# The four-way head-to-head ranking on the common reference scale (RQ4).
relative_efficiency_common_table <- function(out, reference = "M1") {
  relative_efficiency_table(out, reference = reference, scale = NULL,
                            metric = "mse_common", include_m4_forms = FALSE)
}

# Scale gap between VE_h and VE_r per scenario (justifies the common-scale
# comparison; large gaps would warn that the four-way ranking mixes estimands).
scale_gap_summary_table <- function(out) {
  if (is.null(out$config)) return(NULL)
  scale_gap_table(out$config)
}

# M4 misspecification tier: every candidate form scored against every truth,
# signed bias by reporting day (own risk scale).
m4_per_form_bias <- function(out, days = c(30, 60, 90, 120, 150, 180)) {
  pd <- out$per_day[grepl("^M4_", out$per_day$method), ]
  if (nrow(pd) == 0) return(NULL)
  metric_by_day(pd, "signed_bias", days)
}

# Convergence-rate table, wide by method.
convergence_table <- function(out) {
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr is required.")
  c1 <- out$convergence[, c("scenario", "n", "loss_to_followup", "method", "convergence_rate")]
  w <- tidyr::pivot_wider(c1, names_from = "method", values_from = "convergence_rate")
  as.data.frame(w[order(w$scenario, w$n, w$loss_to_followup), ])
}

# M4 LOOIC functional-form selection table: counts and proportions per cell.
m4_selection_table <- function(out) {
  sel <- out$m4_selection
  if (is.null(sel)) return(NULL)
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required.")
  tot <- dplyr::summarise(dplyr::group_by(sel, .data$scenario, .data$n, .data$loss_to_followup),
                          total = sum(.data$Freq), .groups = "drop")
  m <- merge(sel, tot, by = c("scenario", "n", "loss_to_followup"))
  m$proportion <- m$Freq / m$total
  m[order(m$scenario, m$n, m$loss_to_followup, -m$Freq), ]
}

# ---- Rendering --------------------------------------------------------------

# Render a data frame to gt HTML (no headless-browser dependency needed).
render_gt <- function(df, title, subtitle = NULL, file = NULL) {
  if (!requireNamespace("gt", quietly = TRUE)) return(invisible(NULL))
  g <- gt::gt(df)
  g <- gt::tab_header(g, title = title, subtitle = subtitle)
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  if (length(num_cols) > 0) g <- gt::fmt_number(g, columns = num_cols, decimals = 3)
  if (!is.null(file)) gt::gtsave(g, file)
  g
}

# Write the full table set to CSV (+ headline gt HTML). Returns the tables list.
write_summary_tables <- function(out, out_dir = file.path("outputs", "summary"),
                                 days = c(30, 60, 90, 120, 150, 180), render_html = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- list(
    bias = bias_by_day(out$per_day, days),
    rmse = rmse_by_day(out$per_day, days),
    coverage = coverage_by_day(out$per_day, days),
    width = width_by_day(out$per_day, days),
    ise = ise_table(out),
    relative_efficiency = relative_efficiency_table(out),
    relative_efficiency_common = relative_efficiency_common_table(out),
    scale_gap = scale_gap_summary_table(out),
    m4_per_form_bias = m4_per_form_bias(out, days),
    convergence = convergence_table(out),
    m4_selection = m4_selection_table(out)
  )
  for (nm in names(tables)) {
    if (!is.null(tables[[nm]])) {
      write.csv(tables[[nm]], file.path(out_dir, paste0(nm, ".csv")), row.names = FALSE)
    }
  }
  if (render_html) {
    titles <- c(bias = "Signed bias by day", coverage = "95% interval coverage by day",
                convergence = "Convergence rate by method", m4_selection = "M4 LOOIC form selection")
    sub <- "M1-M3/M1b: hazard-scale VE_h; M4: risk-scale interval VE_r. M1b is a supplementary comparator, not one of the four thesis methods."
    for (nm in names(titles)) {
      if (!is.null(tables[[nm]])) {
        try(render_gt(tables[[nm]], titles[[nm]], sub,
                      file = file.path(out_dir, paste0(nm, ".html"))), silent = TRUE)
      }
    }
  }
  invisible(tables)
}
