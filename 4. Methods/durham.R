################################################################################
# DURHAM-STYLE WANING TOOLKIT (Reusable Functions)
# - Single binary treatment (two groups)
# - Outcome: time-to-event with censoring
# - Uses scaled Schoenfeld residual process via cox.zph (identity time scale)
################################################################################

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

#-----------------------------#
# 1) Input validation / prep  #
#-----------------------------#
durham_prep_data <- function(data,
                             time_col   = "time",
                             status_col = "status",
                             treat_col  = "treat",
                             ref_level  = "placebo",
                             vax_level  = "nirsevimab") {
  stopifnot(is.data.frame(data))

  if (!all(c(time_col, status_col, treat_col) %in% names(data))) {
    stop("Data must contain columns: ", paste(c(time_col, status_col, treat_col), collapse = ", "))
  }

  out <- data %>%
    transmute(
      time   = as.numeric(.data[[time_col]]),
      status = as.integer(.data[[status_col]]),
      treat  = as.character(.data[[treat_col]])
    )

  # Ensure binary treatment with explicit reference
  out$treat <- factor(out$treat, levels = c(ref_level, vax_level))

  if (any(is.na(out$treat))) {
    stop("Some rows have treat not equal to ref_level or vax_level. ",
         "Check ref_level/vax_level or your data coding.")
  }

  if (!all(out$status %in% c(0, 1))) {
    stop("status must be coded 0/1.")
  }

  if (any(out$time <= 0, na.rm = TRUE)) {
    warning("Some time values are <= 0. Ensure time is strictly positive for Cox models.")
  }

  out
}

#-----------------------------#
# 2) Fit Cox + basic summary  #
#-----------------------------#
durham_fit_cox <- function(dat, ties = "efron") {
  fit <- coxph(Surv(time, status) ~ treat, data = dat, ties = ties, x = TRUE)

  hr <- unname(exp(coef(fit))[1])
  ve <- (1 - hr) * 100

  list(
    fit = fit,
    hr = hr,
    ve = ve,
    n = nrow(dat),
    events = sum(dat$status)
  )
}

#----------------------------------------------#
# 3) PH / waning test (Grambsch–Therneau test) #
#----------------------------------------------#
durham_ph_test <- function(cox_fit, transform = "identity") {
  zph <- cox.zph(cox_fit, transform = transform)
  # For single-covariate model, first row corresponds to treat
  p_treat <- unname(zph$table[1, "p"])
  p_global <- unname(zph$table[nrow(zph$table), "p"])

  list(
    zph = zph,
    p_treat = p_treat,
    p_global = p_global
  )
}

#------------------------------------------------------------#
# 4) Durham curve: construct beta(t) and VE(t) + smoothing    #
#------------------------------------------------------------#
durham_curve <- function(dat,
                         t_grid = NULL,
                         smooth = c("loess", "ns"),
                         span = 0.75,
                         ns_df = 4,
                         zph_transform = "identity",
                         loess_degree = 1) {
  smooth <- match.arg(smooth)

  # Fit Cox
  fit_obj <- durham_fit_cox(dat)
  fit <- fit_obj$fit
  beta_hat <- unname(coef(fit)[1])

  # Scaled Schoenfeld residual process on chosen time transform
  zph <- cox.zph(fit, transform = zph_transform)

  # Time points for residual process
  resid_df <- tibble(
    time = as.numeric(zph$x),
    beta_t_raw = beta_hat + as.numeric(zph$y[, 1])
  ) %>% arrange(time)

  # Prediction grid
  if (is.null(t_grid)) {
    t0 <- max(min(dat$time, na.rm=TRUE), 1e-6)
    t_grid <- seq(t0, max(dat$time, na.rm=TRUE), length.out = 200)

  }

  # Smooth beta(t)
  if (smooth == "loess") {
    lo_fit <- loess(
      beta_t_raw ~ time,
      data = resid_df,
      span = span,
      degree = loess_degree,
      control = loess.control(surface = "direct")
    )
    beta_s <- predict(lo_fit, newdata = data.frame(time = t_grid))
  } else {
    ns_fit <- lm(beta_t_raw ~ splines::ns(time, df = ns_df), data = resid_df)
    beta_s <- predict(ns_fit, newdata = data.frame(time = t_grid))
  }

  ve_df <- tibble(
    time = t_grid,
    beta_t = beta_s,
    HR_t = exp(beta_s),
    VE_t = (1 - HR_t) * 100
  )

  list(
    cox = fit_obj,
    zph = zph,
    residual_process = resid_df,
    curve = ve_df
  )
}

#------------------------------------------------------------#
# 5) Stratified bootstrap bands for VE(t) (pointwise 95% CI)  #
#------------------------------------------------------------#
durham_bootstrap_bands <- function(dat,
                                   t_grid = NULL,
                                   B = 500,
                                   smooth = c("loess", "ns"),
                                   span = 0.75,
                                   ns_df = 4,
                                   zph_transform = "identity",
                                   seed = 42) {
  smooth <- match.arg(smooth)
  set.seed(seed)

  if (is.null(t_grid)) {
    t_grid <- seq(0, max(dat$time), length.out = 200)
  }

  # Stratified bootstrap (version-proof)
  bootstrap_sample <- function(d) {
    d %>%
      group_by(treat) %>%
      sample_frac(size = 1, replace = TRUE) %>%
      ungroup()
  }

  boot_mat <- matrix(NA_real_, nrow = length(t_grid), ncol = B)

  for (b in seq_len(B)) {
    db <- bootstrap_sample(dat)

    out <- try(
      durham_curve(
        db, t_grid = t_grid,
        smooth = smooth, span = span, ns_df = ns_df,
        zph_transform = zph_transform
      )$curve$VE_t,
      silent = TRUE
    )

    if (!inherits(out, "try-error")) {
      boot_mat[, b] <- out
    }
  }

  # Pointwise 95% bands
  ve_lower <- apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
  ve_upper <- apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE)

  list(
    t_grid = t_grid,
    boot_mat = boot_mat,
    ve_lower = ve_lower,
    ve_upper = ve_upper,
    n_success = sum(colSums(!is.na(boot_mat)) > 0)
  )
}

#--------------------------------#
# 6) One-shot full analysis       #
#--------------------------------#
durham_waning_analysis <- function(data,
                                   time_col = "time",
                                   status_col = "status",
                                   treat_col = "treat",
                                   ref_level = "placebo",
                                   vax_level = "nirsevimab",
                                   t_grid = NULL,
                                   smooth = "loess",
                                   span = 0.75,
                                   ns_df = 4,
                                   zph_transform = "identity",
                                   B = 500,
                                   seed = 42) {
  # Prep
  dat <- durham_prep_data(
    data,
    time_col = time_col, status_col = status_col, treat_col = treat_col,
    ref_level = ref_level, vax_level = vax_level
  )

  # Main curve
  main <- durham_curve(
    dat, t_grid = t_grid,
    smooth = smooth, span = span, ns_df = ns_df,
    zph_transform = zph_transform
  )

  # PH test
  ph <- durham_ph_test(main$cox$fit, transform = zph_transform)

  # Bootstrap bands
  boot <- durham_bootstrap_bands(
    dat, t_grid = main$curve$time,
    B = B, smooth = smooth, span = span, ns_df = ns_df,
    zph_transform = zph_transform,
    seed = seed
  )

  # Merge bands
  curve_ci <- main$curve %>%
    mutate(
      VE_lower = boot$ve_lower,
      VE_upper = boot$ve_upper
    )

  list(
    data = dat,
    cox = main$cox,
    ph = ph,
    residual_process = main$residual_process,
    curve = curve_ci,
    bootstrap = boot,
    settings = list(
      time_col = time_col, status_col = status_col, treat_col = treat_col,
      ref_level = ref_level, vax_level = vax_level,
      smooth = smooth, span = span, ns_df = ns_df,
      zph_transform = zph_transform,
      B = B, seed = seed
    )
  )
}

#-------------------------------#
# 7) Plot helper                #
#-------------------------------#
plot_durham_ve <- function(res, ylim = c(0, 100)) {
  ve0 <- res$cox$ve

  ggplot(res$curve, aes(x = time, y = VE_t)) +
    geom_ribbon(aes(ymin = VE_lower, ymax = VE_upper), alpha = 0.25) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = ve0, linetype = "dashed") +
    coord_cartesian(ylim = ylim) +
    labs(
      title = "Time-varying VE(t) (Durham-style) with bootstrap 95% bands",
      subtitle = paste0("B = ", res$settings$B,
                        "; dashed line = overall Cox VE (", round(ve0, 1), "%)"),
      x = "Days since dose",
      y = "Vaccine efficacy (%)"
    ) +
    theme_bw()
}


#------------------------------------------------------------#
# 4b) Durham curve using SPLINES (explicit spline smoother)   #
#     - Smoothes Schoenfeld residual process with splines     #
#------------------------------------------------------------#
durham_curve_spline <- function(dat,
                                t_grid = NULL,
                                spline = c("ns", "bs"),
                                df = 4,
                                degree = 3,         # for bs()
                                intercept = FALSE,  # for bs()
                                zph_transform = "identity",
                                ties = "efron") {

  spline <- match.arg(spline)

  # Fit Cox (overall effect)
  fit_obj <- durham_fit_cox(dat, ties = ties)
  fit <- fit_obj$fit
  beta_hat <- unname(coef(fit)[1])

  # Scaled Schoenfeld residual process
  zph <- cox.zph(fit, transform = zph_transform)

  resid_df <- tibble::tibble(
    time = as.numeric(zph$x),
    beta_t_raw = beta_hat + as.numeric(zph$y[, 1])
  ) %>%
    dplyr::arrange(time)

  # Prediction grid (avoid 0)
  if (is.null(t_grid)) {
    t0 <- max(min(dat$time, na.rm = TRUE), 1e-6)
    t_grid <- seq(t0, max(dat$time, na.rm = TRUE), length.out = 200)
  }

  # Fit spline smoother
  if (spline == "ns") {
    spline_fit <- lm(beta_t_raw ~ splines::ns(time, df = df), data = resid_df)
  } else {
    spline_fit <- lm(beta_t_raw ~ splines::bs(time, df = df, degree = degree, intercept = intercept),
                     data = resid_df)
  }

  beta_s <- predict(spline_fit, newdata = data.frame(time = t_grid))

  ve_df <- tibble::tibble(
    time = t_grid,
    beta_t = beta_s,
    HR_t = exp(beta_s),
    VE_t = (1 - HR_t) * 100
  )

  list(
    cox = fit_obj,
    zph = zph,
    residual_process = resid_df,
    spline_fit = spline_fit,
    curve = ve_df,
    settings = list(
      spline = spline,
      df = df,
      degree = degree,
      intercept = intercept,
      zph_transform = zph_transform,
      ties = ties
    )
  )
}


#------------------------------------------------------------#
# 5b) Bootstrap bands for SPLINE Durham curve (pointwise CI)  #
#------------------------------------------------------------#
durham_bootstrap_bands_spline <- function(dat,
                                         t_grid = NULL,
                                         B = 500,
                                         spline = c("ns", "bs"),
                                         df = 4,
                                         degree = 3,
                                         intercept = FALSE,
                                         zph_transform = "identity",
                                         ties = "efron",
                                         seed = 42) {

  spline <- match.arg(spline)
  set.seed(seed)

  if (is.null(t_grid)) {
    t0 <- max(min(dat$time, na.rm = TRUE), 1e-6)
    t_grid <- seq(t0, max(dat$time, na.rm = TRUE), length.out = 200)
  }

  bootstrap_sample <- function(d) {
    d %>%
      dplyr::group_by(treat) %>%
      dplyr::sample_frac(size = 1, replace = TRUE) %>%
      dplyr::ungroup()
  }

  boot_mat <- matrix(NA_real_, nrow = length(t_grid), ncol = B)

  for (b in seq_len(B)) {
    db <- bootstrap_sample(dat)

    out <- try(
      durham_curve_spline(
        db, t_grid = t_grid,
        spline = spline, df = df,
        degree = degree, intercept = intercept,
        zph_transform = zph_transform,
        ties = ties
      )$curve$VE_t,
      silent = TRUE
    )

    if (!inherits(out, "try-error")) boot_mat[, b] <- out
  }

  list(
    t_grid = t_grid,
    boot_mat = boot_mat,
    ve_lower = apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    ve_upper = apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE),
    n_success = sum(colSums(!is.na(boot_mat)) > 0)
  )
}



# ==============================================================================
# 2) SIMULTANEOUS BOOTSTRAP BANDS (Durham and/or TDC) ✅
#    Instead of pointwise quantiles, controls the max deviation across time.
#    Band: VE_main(t) ± c where c = 95th percentile of max|VE_b(t)-VE_main(t)|.
#
#    You can use this for:
#     - Durham: pass boot_mat from durham_bootstrap_bands()
#     - TDC:    pass boot_mat from tdc_bootstrap_bands_logt()
#
#    NOTE: This yields symmetric simultaneous bands around the main curve.
# ==============================================================================

bootstrap_simultaneous_bands <- function(main_vec,
                                        boot_mat,
                                        level = 0.95,
                                        na.rm = TRUE) {
  stopifnot(is.numeric(main_vec))
  stopifnot(is.matrix(boot_mat))
  stopifnot(nrow(boot_mat) == length(main_vec))

  # Drop bootstrap columns with too many NA if desired
  good_cols <- rep(TRUE, ncol(boot_mat))
  if (na.rm) {
    good_cols <- colSums(is.na(boot_mat)) < nrow(boot_mat)
  }
  bm <- boot_mat[, good_cols, drop = FALSE]
  if (ncol(bm) == 0) stop("No usable bootstrap replicates after NA filtering.")

  # Deviations from main
  dev <- sweep(bm, 1, main_vec, FUN = "-")

  # Max absolute deviation per replicate
  max_abs_dev <- apply(dev, 2, function(v) max(abs(v), na.rm = TRUE))

  alpha <- 1 - level
  c_val <- unname(quantile(max_abs_dev, probs = 1 - alpha, na.rm = TRUE))

  list(
    band_const = c_val,
    lower = main_vec - c_val,
    upper = main_vec + c_val,
    level = level,
    n_success = ncol(bm)
  )
}

# Durham usage:
# main_vec <- res$curve$VE_t
# sim <- bootstrap_simultaneous_bands(main_vec, res$bootstrap$boot_mat, level=0.95)
# res$curve <- res$curve %>% mutate(VE_lower_sim = sim$lower, VE_upper_sim = sim$upper)

# TDC usage:
# main_vec <- res_tdc$curve$VE_t
# sim <- bootstrap_simultaneous_bands(main_vec, res_tdc$boot$boot_mat, level=0.95)
# res_tdc$curve <- res_tdc$curve %>% mutate(VE_lower_sim = sim$lower, VE_upper_sim = sim$upper)


# ==============================================================================
# 3) m-out-of-n BOOTSTRAP FOR DURHAM ✅ (reduces explosions)
#    Resample m < n (within treatment strata), fit Durham curve each time.
#    Typically m = floor(m_frac * n), e.g. m_frac=0.8.
#
#    NOTE: This is an *approximate* stabilization trick. It often reduces extreme
#    bootstrap replicates for rare-event endpoints.
# ==============================================================================

durham_bootstrap_bands_moutofn <- function(dat,
                                          t_grid = NULL,
                                          B = 500,
                                          m_frac = 0.8,
                                          smooth = c("loess", "ns"),
                                          span = 0.75,
                                          ns_df = 4,
                                          zph_transform = "identity",
                                          loess_degree = 1,
                                          seed = 42) {
  smooth <- match.arg(smooth)
  set.seed(seed)

  if (is.null(t_grid)) {
    t0 <- max(min(dat$time, na.rm = TRUE), 1e-6)
    t_grid <- seq(t0, max(dat$time, na.rm = TRUE), length.out = 200)
  }

  n_total <- nrow(dat)
  m_total <- max(2, floor(m_frac * n_total))

  # Allocate m across strata proportional to stratum sizes
  n_by <- dat %>% count(treat, name = "n_stratum")
  n_by <- n_by %>% mutate(m_stratum = pmax(1, round(m_total * n_stratum / sum(n_stratum))))

  bootstrap_sample_m <- function(d) {
    d %>%
      group_by(treat) %>%
      group_modify(function(.x, .y) {
        m_s <- n_by$m_stratum[n_by$treat == .y$treat]
        # sample_n with replace to mimic bootstrap; size = m_s
        dplyr::sample_n(.x, size = m_s, replace = TRUE)
      }) %>%
      ungroup()
  }

  boot_mat <- matrix(NA_real_, nrow = length(t_grid), ncol = B)

  for (b in seq_len(B)) {
    db <- bootstrap_sample_m(dat)

    out <- try(
      durham_curve(
        db,
        t_grid = t_grid,
        smooth = smooth,
        span = span,
        ns_df = ns_df,
        zph_transform = zph_transform,
        loess_degree = loess_degree
      )$curve$VE_t,
      silent = TRUE
    )
    if (!inherits(out, "try-error")) boot_mat[, b] <- out
  }

  ve_lower <- apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
  ve_upper <- apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE)

  list(
    t_grid = t_grid,
    boot_mat = boot_mat,
    ve_lower = ve_lower,
    ve_upper = ve_upper,
    n_success = sum(colSums(!is.na(boot_mat)) > 0),
    m_frac = m_frac,
    B = B
  )
}

# Usage:
# boot_m <- durham_bootstrap_bands_moutofn(dat, t_grid=curve_main$time, B=500, m_frac=0.8)
# curve_ci <- curve_main %>% mutate(VE_lower_m = boot_m$ve_lower, VE_upper_m = boot_m$ve_upper)
