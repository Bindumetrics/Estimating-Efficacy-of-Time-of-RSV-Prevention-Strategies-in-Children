# ---------------------------------------------------------------------------
# Figures: pointwise metric curves, truth vs estimate, replicate spaghetti.
#
# Pointwise curves are the point of the whole design -- they reveal WHERE a
# method fails, which a grid-averaged scalar hides entirely.
#
# Clipping: VE is NEVER clipped internally. It is clipped only here, for
# plotting, and the clipped fraction is printed so the reader knows how much of
# the tail is being hidden by the axis limits.
#
#   Rscript scripts/figures.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(ggplot2) })
source("R/dgp.R"); source("R/scenarios.R")

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

pw <- read.csv("outputs/pointwise_metrics.csv", stringsAsFactors = FALSE)
raw_path <- if (file.exists("outputs/raw_estimates.parquet"))
  "outputs/raw_estimates.parquet" else "outputs/raw_estimates.rds"
est <- if (grepl("parquet$", raw_path)) arrow::read_parquet(raw_path) else readRDS(raw_path)

METHOD_LEVELS <- c("cox_tdc", "durham", "tian", "pooled_logistic", "bayes")
METHOD_LABELS <- c("Cox-TDC", "Durham", "Tian", "Pooled logistic", "Bayesian")
pw$method <- factor(pw$method, METHOD_LEVELS, METHOD_LABELS)
est$method <- factor(est$method, METHOD_LEVELS, METHOD_LABELS)

# Colour-blind-safe qualitative palette, consistent across every figure so a
# method keeps the same colour throughout the chapter.
PAL <- c("Cox-TDC" = "#4C6EF5", "Durham" = "#E8590C", "Tian" = "#0CA678",
         "Pooled logistic" = "#AE3EC9", "Bayesian" = "#F59F00")

base <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom", legend.title = element_blank())

# --- 1. pointwise metric curves, one panel per metric, faceted by scenario ---
metrics <- list(
  bias     = list(col = "bias",     lab = "Bias(t)",            ref = 0),
  rmse     = list(col = "rmse",     lab = "RMSE(t)",            ref = NA),
  coverage = list(col = "coverage", lab = "Coverage(t)",        ref = 0.95),
  width    = list(col = "width",    lab = "Mean CI width(t)",   ref = NA)
)

for (mn in names(metrics)) {
  m <- metrics[[mn]]
  d <- pw[!is.na(pw[[m$col]]), ]
  mcse_col <- paste0(m$col, "_mcse")
  p <- ggplot(d, aes(t, .data[[m$col]], colour = method, fill = method)) +
    { if (mcse_col %in% names(d))
        geom_ribbon(aes(ymin = .data[[m$col]] - 1.96 * .data[[mcse_col]],
                        ymax = .data[[m$col]] + 1.96 * .data[[mcse_col]]),
                    colour = NA, alpha = 0.15) } +
    { if (is.finite(m$ref))
        geom_hline(yintercept = m$ref, linetype = "22", colour = "grey40") } +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ scenario, scales = "free_y") +
    scale_colour_manual(values = PAL) + scale_fill_manual(values = PAL) +
    labs(x = "Days since randomisation", y = m$lab,
         subtitle = paste0(m$lab,
           " by method and scenario. Bands are 1.96 x Monte Carlo SE.")) +
    base
  ggsave(sprintf("outputs/figures/pointwise_%s.png", mn), p,
         width = 11, height = 8, dpi = 150)
  say("wrote outputs/figures/pointwise_%s.png", mn)
}

# --- 2. truth vs mean estimate ----------------------------------------------
tv <- pw[!is.na(pw$bias), ]
tv$est_mean <- tv$ve_true + tv$bias
p <- ggplot(tv, aes(t)) +
  geom_line(aes(y = ve_true), colour = "grey20", linewidth = 0.9) +
  geom_line(aes(y = est_mean, colour = method), linewidth = 0.6) +
  facet_wrap(~ scenario) +
  scale_colour_manual(values = PAL) +
  labs(x = "Days since randomisation", y = expression(VE[h](t)),
       subtitle = "Grey = true VE_h(t) (the estimand); coloured = mean estimate across replicates") +
  base
ggsave("outputs/figures/truth_vs_estimate.png", p, width = 11, height = 8, dpi = 150)
say("wrote outputs/figures/truth_vs_estimate.png")

# --- 3. replicate spaghetti, reference scenario -----------------------------
YLIM <- c(-0.5, 1)
sp <- est[est$scenario == "ref" & est$rep <= 100 & !is.na(est$ve_hat), ]
if (nrow(sp)) {
  clipped <- mean(sp$ve_hat < YLIM[1] | sp$ve_hat > YLIM[2])
  say("spaghetti plot: %.2f%% of estimates fall outside the [%.1f, %.1f] axis and",
      100 * clipped, YLIM[1], YLIM[2])
  say("  are clipped FOR DISPLAY ONLY -- no internal clipping is applied anywhere.")
  tru <- pw[pw$scenario == "ref", c("t", "ve_true", "method")]
  p <- ggplot(sp, aes(t, ve_hat, group = rep)) +
    geom_line(alpha = 0.08, colour = "#4C6EF5") +
    geom_line(data = unique(tru[, c("t", "ve_true")]),
              aes(t, ve_true, group = 1), colour = "black", linewidth = 0.9) +
    facet_wrap(~ method, nrow = 1) +
    coord_cartesian(ylim = YLIM) +
    labs(x = "Days since randomisation", y = expression(hat(VE)(t)),
         subtitle = sprintf(
           "First 100 replicates, reference scenario. Black = truth. %.1f%% of estimates clipped for display.",
           100 * clipped)) +
    base
  ggsave("outputs/figures/spaghetti_ref.png", p, width = 12, height = 3.6, dpi = 150)
  say("wrote outputs/figures/spaghetti_ref.png")
}

# --- 4. coverage vs width trade-off -----------------------------------------
sm <- read.csv("outputs/summary_metrics.csv", stringsAsFactors = FALSE)
sm$method <- factor(sm$method, METHOD_LEVELS, METHOD_LABELS)
p <- ggplot(sm, aes(mean_width, overall_coverage, colour = method)) +
  geom_hline(yintercept = 0.95, linetype = "22", colour = "grey40") +
  geom_point(size = 2) +
  scale_colour_manual(values = PAL) +
  labs(x = "Mean interval width", y = "Overall coverage",
       subtitle = "A method below the dashed line buys its narrow intervals with under-coverage") +
  base
ggsave("outputs/figures/coverage_vs_width.png", p, width = 7, height = 5, dpi = 150)
say("wrote outputs/figures/coverage_vs_width.png")
