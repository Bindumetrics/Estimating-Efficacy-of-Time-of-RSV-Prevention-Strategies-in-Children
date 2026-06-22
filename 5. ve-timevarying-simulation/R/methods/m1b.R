# methods/m1b.R
# -----------------------------------------------------------------------------
# M1b: flexible Royston-Parmar / spline comparator using rstpm2.
#
# This is a simulation-only supplementary comparator, not one of the four thesis
# methods. It is included so M1 is not the only restrictive Cox-style method
# being compared with flexible alternatives.
#
# Model form:
#   rstpm2::stpm2(Surv(time, status) ~ z, df = baseline_df, tvc = list(z = tvc_df))
#
# The treatment effect is allowed to vary smoothly with log time. VE_h(t) is
# computed from the fitted hazard ratio h_v(t) / h_c(t) on the common day grid.
# -----------------------------------------------------------------------------

fit_m1b <- function(data, grid, baseline_df = 4L, tvc_df = 3L, conf_level = 0.95) {
  if (!requireNamespace("rstpm2", quietly = TRUE)) {
    return(m1b_failure("The rstpm2 package is required for M1b."))
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    return(m1b_failure("The survival package is required for M1b."))
  }
  suppressPackageStartupMessages(require("rstpm2", quietly = TRUE, character.only = TRUE))

  required <- c("time", "status", "z")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    return(m1b_failure(paste("Missing required columns:", paste(missing, collapse = ", "))))
  }
  if (sum(data$status == 1L, na.rm = TRUE) < baseline_df + tvc_df + 2L) {
    return(m1b_failure("Too few events for the requested spline degrees of freedom."))
  }

  data <- data[, required]
  data$z <- as.numeric(data$z)

  fit <- tryCatch(
    rstpm2::stpm2(
      survival::Surv(time, status) ~ z,
      data = data,
      df = as.integer(baseline_df),
      tvc = list(z = as.integer(tvc_df))
    ),
    warning = function(w) structure(list(warning = conditionMessage(w)), class = "m1b_warning"),
    error = function(e) structure(list(error = conditionMessage(e)), class = "m1b_error")
  )

  if (inherits(fit, "m1b_warning")) {
    warning_message <- fit$warning
    fit2 <- tryCatch(
      suppressWarnings(rstpm2::stpm2(
        survival::Surv(time, status) ~ z,
        data = data,
        df = as.integer(baseline_df),
        tvc = list(z = as.integer(tvc_df))
      )),
      error = function(e) structure(list(error = e$message), class = "m1b_error")
    )
    if (inherits(fit2, "m1b_error")) return(m1b_failure(fit2$error))
    fit <- fit2
  } else if (inherits(fit, "m1b_error")) {
    return(m1b_failure(fit$error))
  } else {
    warning_message <- NA_character_
  }

  estimates <- predict_m1b_rstpm2(fit, grid = grid, conf_level = conf_level)
  if (is.null(estimates) || any(!is.finite(estimates$ve))) {
    return(m1b_failure("M1b prediction failed or returned non-finite estimates."))
  }

  list(
    method = "M1b_supplementary",
    label = "M1b, supplementary rstpm2 flexible parametric comparator",
    converged = TRUE,
    failure_reason = NA_character_,
    warning = warning_message,
    fit = fit,
    baseline_df = as.integer(baseline_df),
    tvc_df = as.integer(tvc_df),
    estimates = estimates
  )
}

predict_m1b <- function(fit_object, grid, conf_level = 0.95) {
  if (!isTRUE(fit_object$converged)) return(fit_object)
  predict_m1b_rstpm2(fit_object$fit, grid = grid, conf_level = conf_level)
}

predict_m1b_rstpm2 <- function(fit, grid, conf_level = 0.95) {
  stopifnot(all(grid > 0))
  alpha <- 1 - conf_level

  nd_vaccine <- data.frame(time = grid, z = 1)
  nd_placebo <- data.frame(time = grid, z = 0)

  haz_v <- rstpm2::predict(fit, newdata = nd_vaccine, type = "hazard", se.fit = TRUE, level = conf_level)
  haz_c <- rstpm2::predict(fit, newdata = nd_placebo, type = "hazard", se.fit = TRUE, level = conf_level)

  haz_v <- as.data.frame(haz_v)
  haz_c <- as.data.frame(haz_c)
  names(haz_v) <- tolower(names(haz_v))
  names(haz_c) <- tolower(names(haz_c))

  required <- c("estimate", "lower", "upper")
  if (!all(required %in% names(haz_v)) || !all(required %in% names(haz_c))) {
    stop("rstpm2 hazard predictions did not include estimate/lower/upper columns.")
  }

  eps <- .Machine$double.eps
  hr <- pmax(haz_v$estimate, eps) / pmax(haz_c$estimate, eps)

  # The supplementary comparator uses paired hazard predictions. rstpm2 provides
  # pointwise intervals for each hazard; combining them on the log scale gives a
  # conservative pointwise interval for the hazard ratio without changing the
  # core four-method comparison.
  se_log_hv <- interval_to_log_se(haz_v$lower, haz_v$upper, conf_level)
  se_log_hc <- interval_to_log_se(haz_c$lower, haz_c$upper, conf_level)
  se_log_hr <- sqrt(se_log_hv^2 + se_log_hc^2)
  zcrit <- stats::qnorm(1 - alpha / 2)

  log_hr <- log(hr)
  lower_log_hr <- log_hr - zcrit * se_log_hr
  upper_log_hr <- log_hr + zcrit * se_log_hr

  data.frame(
    method = "M1b_supplementary",
    day = grid,
    log_hr = log_hr,
    se_log_hr = se_log_hr,
    hazard_vaccine = haz_v$estimate,
    hazard_placebo = haz_c$estimate,
    ve = 1 - exp(log_hr),
    lower = 1 - exp(upper_log_hr),
    upper = 1 - exp(lower_log_hr),
    conf_level = conf_level,
    supplementary = TRUE
  )
}

interval_to_log_se <- function(lower, upper, conf_level) {
  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
  (log(pmax(upper, .Machine$double.eps)) - log(pmax(lower, .Machine$double.eps))) / (2 * zcrit)
}

m1b_failure <- function(reason) {
  list(
    method = "M1b_supplementary",
    label = "M1b, supplementary rstpm2 flexible parametric comparator",
    converged = FALSE,
    failure_reason = reason,
    warning = NA_character_,
    fit = NULL,
    baseline_df = NA_integer_,
    tvc_df = NA_integer_,
    estimates = NULL
  )
}

