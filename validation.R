# ======================================================================
# VALIDATION SUITE FOR RECONSTRUCTED IPD BUNDLES
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(glue)
  library(ggplot2)
  library(survival)
  library(readr)
})

# -----------------------------
# Helpers
# -----------------------------
rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))
mae  <- function(x, y) mean(abs(x - y), na.rm = TRUE)
maxae <- function(x, y) max(abs(x - y), na.rm = TRUE)

# Return survival S(t) from survfit at arbitrary times (piecewise-constant, right-continuous)
sf_surv_at <- function(sf, times) {
  s <- summary(sf, times = times, extend = TRUE)
  tibble(time = s$time, surv = s$surv, strata = s$strata)
}

# Clean strata label
clean_strata <- function(x, prefix = "treat=") sub(paste0("^", prefix), "", x)

# -----------------------------
# Core: KM pointwise metrics
# -----------------------------
km_pointwise_metrics <- function(bundle) {

  ipd <- bundle$ipd
  km_clean <- bundle$km_clean

  # survfit from reconstructed IPD
  sf <- survfit(Surv(time, status) ~ treat, data = ipd)

  # Digitised KM points (cleaned) as the reference
  dig <- bind_rows(
    km_clean$placebo    %>% mutate(treat = "placebo"),
    km_clean$nirsevimab %>% mutate(treat = "nirsevimab")
  ) %>%
    arrange(treat, time)

  # Evaluate reconstructed KM at digitised time points
  sf_at <- sf_surv_at(sf, dig$time) %>%
    mutate(treat = clean_strata(strata, "treat=")) %>%
    select(time, treat, surv_hat = surv)

  comp <- dig %>%
    left_join(sf_at, by = c("time", "treat")) %>%
    mutate(err = surv_hat - surv)

  # Metrics per arm
  metrics <- comp %>%
    group_by(endpoint = bundle$endpoint, treat) %>%
    summarise(
      n_points = n(),
      RMSE = rmse(surv_hat, surv),
      MAE  = mae(surv_hat, surv),
      MaxAE = maxae(surv_hat, surv),
      .groups = "drop"
    )

  # KS test on survival probabilities at digitised points (per arm)
  # (Note: This is a pragmatic check; you used a KS test already from IPDfromKM;
  # here we test digitised vs reconstructed-at-digitised-times.)
  ks <- comp %>%
    group_by(endpoint = bundle$endpoint, treat) %>%
    summarise(
      D = suppressWarnings(ks.test(surv_hat, surv)$statistic[[1]]),
      p_value = suppressWarnings(ks.test(surv_hat, surv)$p.value),
      .groups = "drop"
    )

  list(comp = comp, metrics = metrics, ks = ks, survfit = sf)
}

# -----------------------------
# Numbers-at-risk validation
# -----------------------------
risk_table_validation <- function(bundle, survfit_obj = NULL) {

  ipd <- bundle$ipd
  trisk <- bundle$trisk
  nrisk_reported <- bundle$nrisk

  sf <- if (is.null(survfit_obj)) survfit(Surv(time, status) ~ treat, data = ipd) else survfit_obj

  rs <- summary(sf, times = trisk, extend = TRUE)

  risk_hat <- tibble(
    endpoint = bundle$endpoint,
    time = rs$time,
    treat = clean_strata(rs$strata, "treat="),
    n_risk_hat = rs$n.risk
  )

  risk_rep <- tibble(
    endpoint = bundle$endpoint,
    time = trisk
  ) %>%
    mutate(
      placebo_reported = nrisk_reported$placebo,
      nirsevimab_reported = nrisk_reported$nirsevimab
    ) %>%
    pivot_longer(cols = c(placebo_reported, nirsevimab_reported),
                 names_to = "treat",
                 values_to = "n_risk_reported") %>%
    mutate(
      treat = case_when(
        treat == "placebo_reported" ~ "placebo",
        treat == "nirsevimab_reported" ~ "nirsevimab",
        TRUE ~ treat
      )
    )

  comp <- risk_rep %>%
    left_join(risk_hat, by = c("endpoint", "time", "treat")) %>%
    mutate(diff = n_risk_hat - n_risk_reported)

  summary <- comp %>%
    group_by(endpoint, treat) %>%
    summarise(
      max_abs_diff = max(abs(diff), na.rm = TRUE),
      any_mismatch = any(diff != 0, na.rm = TRUE),
      .groups = "drop"
    )

  list(comp = comp, summary = summary)
}

# -----------------------------
# Event count validation (optional)
# -----------------------------
# Provide published cumulative event counts if you have them.
# Format:
# published_counts <- list(
#   MD_RSV_LRTI = tibble(time = c(150), placebo = c(51), nirsevimab = c(19)),
#   RSV_Hospitalization = tibble(time = c(150), placebo = c(21), nirsevimab = c(9)),
#   Severe_RSV_LRTI = tibble(time = c(90, 180), placebo = c(33, 62), nirsevimab = c(6, 19)),
#   Any_RSV_LRTI = tibble(time = c(90, 180), placebo = c(...), nirsevimab = c(...))
# )

event_count_validation <- function(bundle, published_counts = NULL, default_times = NULL) {

  ipd <- bundle$ipd
  endpoint <- bundle$endpoint

  # Decide timepoints:
  if (!is.null(published_counts) && endpoint %in% names(published_counts)) {
    pub <- published_counts[[endpoint]]
    times <- pub$time
  } else {
    # fallback: validate at trisk or user-supplied default_times
    times <- if (!is.null(default_times)) default_times else bundle$trisk
    pub <- NULL
  }

  rec <- map_dfr(times, function(tt) {
    ipd %>%
      mutate(treat = as.character(treat)) %>%
      group_by(treat) %>%
      summarise(
        cum_events = sum(status == 1 & time <= tt, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(time = tt, endpoint = endpoint)
  }) %>%
    select(endpoint, time, treat, cum_events) %>%
    pivot_wider(names_from = treat, values_from = cum_events)

  if (is.null(pub)) {
    return(list(reconstructed = rec, comparison = NULL))
  }

  comp <- pub %>%
    mutate(endpoint = endpoint) %>%
    left_join(rec, by = c("endpoint", "time"), suffix = c("_published", "_reconstructed")) %>%
    mutate(
      diff_placebo = placebo_reconstructed - placebo_published,
      diff_nirsevimab = nirsevimab_reconstructed - nirsevimab_published
    )

  list(reconstructed = rec, comparison = comp)
}

# -----------------------------
# Cox model validation
# -----------------------------
cox_validation <- function(bundle) {
  ipd <- bundle$ipd
  endpoint <- bundle$endpoint

  fit <- coxph(Surv(time, status) ~ treat, data = ipd)
  s <- summary(fit)

  hr <- unname(exp(coef(fit)))[1]
  ci <- exp(confint(fit))[1, ]

  tibble(
    endpoint = endpoint,
    n = nrow(ipd),
    events = sum(ipd$status == 1),
    HR = hr,
    CI_low = ci[1],
    CI_high = ci[2],
    p_value = s$coefficients[1, "Pr(>|z|)"]
  )
}

# -----------------------------
# Overlay plot (per endpoint)
# -----------------------------
plot_overlay_km <- function(bundle, km_metrics_obj = NULL) {
  endpoint <- bundle$endpoint

  km_out <- if (is.null(km_metrics_obj)) km_pointwise_metrics(bundle) else km_metrics_obj
  sf <- km_out$survfit

  # Step curve data
  sf_df <- data.frame(
    time = sf$time,
    surv = sf$surv,
    treat = rep(names(sf$strata), sf$strata)
  ) %>%
    mutate(treat = clean_strata(treat, "treat="))

  dig <- bind_rows(
    bundle$km_clean$placebo    %>% mutate(treat = "placebo"),
    bundle$km_clean$nirsevimab %>% mutate(treat = "nirsevimab")
  )

  ggplot() +
    geom_step(data = sf_df, aes(time, surv, colour = treat), linewidth = 0.9) +
    geom_point(data = dig, aes(time, surv, colour = treat), alpha = 0.35, size = 1) +
    facet_wrap(~treat) +
    theme_bw() +
    labs(
      title = glue("Validation overlay: reconstructed KM vs digitised KM ({endpoint})"),
      x = "Time",
      y = "Survival probability",
      colour = "Arm"
    )
}

# ======================================================================
# RUN VALIDATION ON ALL SAVED ENDPOINTS
# ======================================================================

run_validation_all <- function(rds_dir,
                               out_dir = file.path(rds_dir, "validation_outputs"),
                               published_counts = NULL,
                               default_event_times = NULL,
                               save_plots = TRUE) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
  stopifnot(length(rds_files) > 0)

  all_km_metrics <- list()
  all_km_ks <- list()
  all_risk_summ <- list()
  all_risk_comp <- list()
  all_cox <- list()
  all_event_comp <- list()
  all_event_rec <- list()

  for (f in rds_files) {
    bundle <- readRDS(f)
    endpoint <- bundle$endpoint

    message(glue("---- Validating: {endpoint} ----"))

    km_out <- km_pointwise_metrics(bundle)
    risk_out <- risk_table_validation(bundle, km_out$survfit)
    cox_out <- cox_validation(bundle)
    evt_out <- event_count_validation(bundle,
                                      published_counts = published_counts,
                                      default_times = default_event_times)

    all_km_metrics[[endpoint]] <- km_out$metrics
    all_km_ks[[endpoint]] <- km_out$ks
    all_risk_summ[[endpoint]] <- risk_out$summary
    all_risk_comp[[endpoint]] <- risk_out$comp
    all_cox[[endpoint]] <- cox_out
    all_event_rec[[endpoint]] <- evt_out$reconstructed
    if (!is.null(evt_out$comparison)) all_event_comp[[endpoint]] <- evt_out$comparison

    if (isTRUE(save_plots)) {
      p <- plot_overlay_km(bundle, km_out)
      ggsave(
        filename = file.path(out_dir, paste0(endpoint, "_KM_overlay.png")),
        plot = p, width = 9, height = 4.8, dpi = 300
      )
    }

    # Save pointwise comparison table for traceability
    readr::write_csv(
      km_out$comp %>% mutate(source_rds = basename(f)),
      file.path(out_dir, paste0(endpoint, "_km_pointwise_comparison.csv"))
    )

    # Save risk comparison table
    readr::write_csv(
      risk_out$comp %>% mutate(source_rds = basename(f)),
      file.path(out_dir, paste0(endpoint, "_risk_comparison.csv"))
    )
  }

  # Bind outputs
  km_metrics_tbl <- bind_rows(all_km_metrics)
  km_ks_tbl      <- bind_rows(all_km_ks)
  risk_summary_tbl <- bind_rows(all_risk_summ)
  cox_tbl        <- bind_rows(all_cox)
  event_rec_tbl  <- bind_rows(all_event_rec)
  event_comp_tbl <- if (length(all_event_comp) > 0) bind_rows(all_event_comp) else NULL

  # Export master summaries
  write_csv(km_metrics_tbl, file.path(out_dir, "summary_km_error_metrics.csv"))
  write_csv(km_ks_tbl,      file.path(out_dir, "summary_km_ks_tests.csv"))
  write_csv(risk_summary_tbl, file.path(out_dir, "summary_risk_table_checks.csv"))
  write_csv(cox_tbl,        file.path(out_dir, "summary_cox_hr.csv"))
  write_csv(event_rec_tbl,  file.path(out_dir, "summary_cumulative_events_reconstructed.csv"))
  if (!is.null(event_comp_tbl)) {
    write_csv(event_comp_tbl, file.path(out_dir, "summary_cumulative_events_comparison.csv"))
  }

  list(
    km_error_metrics = km_metrics_tbl,
    km_ks_tests = km_ks_tbl,
    risk_summary = risk_summary_tbl,
    cox_summary = cox_tbl,
    event_reconstructed = event_rec_tbl,
    event_comparison = event_comp_tbl,
    outputs_dir = out_dir
  )
}


# ------------------------------------------------------------------
# PLOTS FOR KAMPMANN
# ------------------------------------------------------------------
kampann_plot <- function(ipd,
                           time_var = "time",
                           status_var = "status",
                           arm_var = "treat",
                           arm_levels = c("placebo", "RSVpreF"),
                           max_day = 180,
                           breaks = c(0, 30, 60, 90, 120, 150, 180),
                           title = "") {

  dat <- ipd %>%
    transmute(
      time   = .data[[time_var]],
      status = .data[[status_var]],
      arm    = .data[[arm_var]]
    ) %>%
    filter(!is.na(time), !is.na(status), !is.na(arm)) %>%
    mutate(
      time = pmin(time, max_day),
      status = ifelse(time >= max_day, 0, status), # ensure admin censoring at max_day
      arm = factor(arm, levels = arm_levels)
    )

  fit <- survfit(Surv(time, status) ~ arm, data = dat)

  # Build step data: cumulative incidence = 1 - survival
  s <- summary(fit, times = sort(unique(c(0, fit$time[fit$n.event > 0], breaks, max_day))))
  plot_df <- tibble(
    time = s$time,
    surv = s$surv,
    arm  = s$strata
  ) %>%
    mutate(
      arm = sub("^arm=", "", arm),
      cuminc_pct = 100 * (1 - surv)
    )

  # Number at risk at chosen breakpoints
  risk_s <- summary(fit, times = breaks)
  risk_df <- tibble(
    time = risk_s$time,
    n_risk = risk_s$n.risk,
    arm = sub("^arm=", "", risk_s$strata)
  )

  p <- ggplot(plot_df, aes(time, cuminc_pct, color = arm)) +
    geom_step(linewidth = 1) +
    scale_x_continuous(breaks = breaks, limits = c(0, max_day)) +
    labs(
      x = "Days after Birth",
      y = "Cumulative Incidence (%)",
      title = title,
      color = NULL
    ) +
    theme_minimal(base_size = 12)

  list(plot = p, fit = fit, risk = risk_df)
}


# ---------------------------------------------------------
#  Plots for Simoes
# --------------------------------------------------------


simoes_plot <- function(s_data){

ipd <- s_data

# 1) Fit KM and Cox
fit_km  <- survfit(Surv(time, status) ~ treat, data = ipd)
fit_cox <- coxph(Surv(time, status) ~ treat, data = ipd)
s <- summary(fit_cox)

hr_txt <- sprintf(
  "HR %.3f (95%% CI %.3f–%.3f), p=%s",
  exp(coef(fit_cox)),
  exp(confint(fit_cox))[1],
  exp(confint(fit_cox))[2],
  format.pval(s$wald["pvalue"], digits = 2, eps = .0001)
)



# 2) Build KM plot (Simões style)
p <- ggsurvplot(
  fit_km,
  data = ipd,
  censor = TRUE,
  censor.shape = 3, censor.size = 2,
  palette = c("Placebo" = "#6EC1E4", "Nirsevimab" = "#D62828"),
  xlim = c(0, 150),
  break.time.by = 30,
  ylim = c(0.90, 1.00),
  xlab = "Days since dose",
  ylab = "Proportion of infants without event",
  legend.title = "",
  legend.labs = c("Placebo", "Nirsevimab"),
  ggtheme = theme_minimal(base_size = 12)
)

p$plot +
  annotate("text", x = 95, y = 0.905, label = hr_txt, hjust = 0, size = 3.5) -> plt
return(plt)
  
}



