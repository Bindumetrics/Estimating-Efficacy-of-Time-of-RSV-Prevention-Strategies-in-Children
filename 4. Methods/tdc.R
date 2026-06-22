
library(survival)
library(dplyr)
library(tibble)

tdc_fit_logt <- function(data,
                         time_col="time", status_col="status", treat_col="treat",
                         ref_level="placebo", vax_level="nirsevimab",
                         eps = 1e-3, ties="efron") {

  dat <- dplyr::transmute(
    data,
    time   = as.numeric(.data[[time_col]]),
    status = as.integer(.data[[status_col]]),
    treat  = factor(as.character(.data[[treat_col]]),
                    levels = c(ref_level, vax_level))
  )

  if (any(is.na(dat$treat))) stop("Treatment levels not matching ref_level/vax_level.")
  if (!all(dat$status %in% c(0,1))) stop("status must be 0/1.")

  dat <- dplyr::mutate(dat, vax = as.integer(treat == vax_level))

  fit <- coxph(
    Surv(time, status) ~ vax + tt(vax),
    data = dat,
    ties = ties,
    tt = function(x, t, ...) x * log(pmax(t, eps))
  )

  cf <- coef(fit)
  b0 <- unname(cf["vax"])
  b1 <- unname(cf[grep("^tt\\(vax\\)", names(cf))])

  list(
    data = dat,
    fit = fit,
    coef = c(b0 = b0, b1 = b1),
    eps = eps,
    ref_level = ref_level,
    vax_level = vax_level
  )
}


tdc_predict_logt <- function(fit_obj, t_grid) {
  b0 <- fit_obj$coef["b0"]
  b1 <- fit_obj$coef["b1"]
  eps <- fit_obj$eps

  beta_t <- b0 + b1 * log(pmax(t_grid, eps))

  tibble::tibble(
    time = t_grid,
    beta_t = beta_t,
    HR_t = exp(beta_t),
    VE_t = (1 - exp(beta_t)) * 100
  )
}



tdc_bootstrap_bands_logt <- function(data,
                                     time_col="time", status_col="status", treat_col="treat",
                                     ref_level="placebo", vax_level="nirsevimab",
                                     t_grid = NULL,
                                     B = 500, seed = 42,
                                     eps = 1e-3, ties="efron") {

  set.seed(seed)

  # prep once to standard names
  dat0 <- dplyr::transmute(
    data,
    time   = as.numeric(.data[[time_col]]),
    status = as.integer(.data[[status_col]]),
    treat  = factor(as.character(.data[[treat_col]]), levels = c(ref_level, vax_level))
  ) %>% dplyr::mutate(vax = as.integer(treat == vax_level))

  if (is.null(t_grid)) t_grid <- seq(0, max(dat0$time), length.out = 200)

  # version-proof stratified bootstrap
  bootstrap_sample <- function(dat) {
    dat %>%
      dplyr::group_by(treat) %>%
      dplyr::sample_frac(size = 1, replace = TRUE) %>%
      dplyr::ungroup()
  }

  boot_mat <- matrix(NA_real_, nrow = length(t_grid), ncol = B)

  for (b in seq_len(B)) {
    db <- bootstrap_sample(dat0)

    # Fit + predict; skip failures safely
    out <- try({
      fit_b <- coxph(
        Surv(time, status) ~ vax + tt(vax),
        data = db,
        ties = ties,
        tt = function(x, t, ...) x * log(pmax(t, eps))
      )
      cf <- coef(fit_b)
      b0 <- unname(cf["vax"])
      b1 <- unname(cf[grep("^tt\\(vax\\)", names(cf))])
      beta_t <- b0 + b1 * log(pmax(t_grid, eps))
      (1 - exp(beta_t)) * 100
    }, silent = TRUE)

    if (!inherits(out, "try-error")) boot_mat[, b] <- out
  }

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

tdc_logt_analysis <- function(data,
                              time_col="time", status_col="status", treat_col="treat",
                              ref_level="placebo", vax_level="nirsevimab",
                              B = 500, seed = 42,
                              eps = 1e-3) {

  fit_obj <- tdc_fit_logt(
    data, time_col, status_col, treat_col,
    ref_level, vax_level,
    eps = eps
  )

  t_grid <- seq(0, max(fit_obj$data$time), length.out = 200)
  pred <- tdc_predict_logt(fit_obj, t_grid)

  boot <- tdc_bootstrap_bands_logt(
    data, time_col, status_col, treat_col,
    ref_level, vax_level,
    t_grid = t_grid,
    B = B, seed = seed,
    eps = eps
  )

  curve <- pred %>%
    dplyr::mutate(VE_lower = boot$ve_lower,
                  VE_upper = boot$ve_upper)

  list(fit_obj = fit_obj, curve = curve, boot = boot)
}

plot_tdc_logt <- function(res, ylim = c(0, 100)) {
  b0 <- res$fit_obj$coef["b0"]
  b1 <- res$fit_obj$coef["b1"]

  ggplot2::ggplot(res$curve, ggplot2::aes(x = time, y = VE_t)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = VE_lower, ymax = VE_upper), alpha = 0.25) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(
      title = "TDC model: VE(t) with log(time) interaction",
      subtitle = paste0("log HR(t) = b0 + b1*log(t);  b0=",
                        round(b0, 3), ", b1=", round(b1, 3),
                        "; bootstrap bands"),
      x = "Days since dose",
      y = "Vaccine efficacy (%)"
    ) +
    ggplot2::theme_bw()
}
