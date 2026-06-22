# plots.R
# -----------------------------------------------------------------------------
# Stage 5 (figures): per-scenario plots of the true VE(t) against the mean
# estimated curve for each method, with mean interval/credible bands.
#
# Each method is plotted against the truth on ITS OWN scale (M1-M3/M1b vs
# hazard-scale VE_h; M4 vs the interval-conditional risk-scale VE_r) because the
# aggregated `truth` column already carries the scale-correct target per method.
# M1b is labelled as a supplementary comparator, not one of the four thesis
# methods.
# -----------------------------------------------------------------------------

m4_label_methods <- function(method) {
  lev <- c("M1", "M1b_supplementary", "M2", "M3", "M4")
  lab <- c("M1 (Cox x log-time)", "M1b (spline, supplementary)",
           "M2 (smoothed Schoenfeld)", "M3 (local PL)", "M4 (Bayesian binomial)")
  factor(lab[match(method, lev)], levels = lab)
}

# Small-multiples figure for one (scenario, n, loss): truth + mean estimate +
# mean band, faceted by method.
plot_scenario_panels <- function(per_day, scenario, n, loss = 0,
                                 ylim = c(0, 1), config = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  d <- per_day[per_day$scenario == scenario & per_day$n == n &
                 per_day$loss_to_followup == loss, ]
  if (nrow(d) == 0) stop("No rows for scenario ", scenario, ", n ", n, ", loss ", loss)
  d$panel <- m4_label_methods(d$method)

  label <- if (!is.null(config)) config$scenarios[[scenario]]$label else scenario
  ggplot2::ggplot(d, ggplot2::aes(x = .data$day)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$mean_lower, ymax = .data$mean_upper),
                         fill = "#9ecae1", alpha = 0.45) +
    ggplot2::geom_line(ggplot2::aes(y = .data$truth), colour = "black", linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(y = .data$mean_ve), colour = "#08519c", linewidth = 0.9) +
    ggplot2::facet_wrap(~panel) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(
      title = sprintf("Scenario %s: %s  (n = %d, loss = %.0f%%)", scenario, label, n, 100 * loss),
      subtitle = "Black = true VE(t) on the method's scale; blue = mean estimate; band = mean 95% interval. M1b is supplementary.",
      x = "Day", y = "VE(t)"
    ) +
    ggplot2::theme_bw(base_size = 11)
}

# Save one panel figure per scenario at a chosen sample size.
save_scenario_plots <- function(out, out_dir = file.path("outputs", "figures"),
                                n = NULL, loss = 0, width = 9, height = 6, dpi = 160) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  per_day <- out$per_day
  config <- out$config
  scenarios <- unique(per_day$scenario)
  if (is.null(n)) n <- max(per_day$n)
  saved <- character(0)
  for (sc in scenarios) {
    avail_n <- unique(per_day$n[per_day$scenario == sc & per_day$loss_to_followup == loss])
    use_n <- if (n %in% avail_n) n else max(avail_n)
    p <- tryCatch(plot_scenario_panels(per_day, sc, use_n, loss, config = config),
                  error = function(e) NULL)
    if (is.null(p)) next
    f <- file.path(out_dir, sprintf("scenario_%s_n%d.png", sc, use_n))
    ggplot2::ggsave(f, p, width = width, height = height, dpi = dpi)
    saved <- c(saved, f)
  }
  saved
}
