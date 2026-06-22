# methods/m2.R
# -----------------------------------------------------------------------------
# M2: Grambsch-Therneau / Durham smoothed Schoenfeld residual method.
#
# Fit an ordinary Cox model, extract the scaled Schoenfeld residual-based time
# varying coefficient estimates for treatment, smooth them over event time, and
# transform beta(t) to hazard-scale VE_h(t) = 1 - exp(beta(t)).
#
# Intervals are pointwise smoother intervals on the log hazard-ratio scale,
# matching the thesis evaluation requirement for like-for-like coverage.
# -----------------------------------------------------------------------------

fit_m2 <- function(data, grid, basis_dimension = 8L, conf_level = 0.95) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    return(m2_failure("The survival package is required for M2."))
  }
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(m2_failure("The mgcv package is required for M2 smoothing."))
  }
  suppressPackageStartupMessages(require("mgcv", quietly = TRUE, character.only = TRUE))

  required <- c("time", "status", "z")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    return(m2_failure(paste("Missing required columns:", paste(missing, collapse = ", "))))
  }

  n_events <- sum(data$status == 1L, na.rm = TRUE)
  if (n_events < 6L) {
    return(m2_failure("Fewer than six events; Schoenfeld smoothing is not stable."))
  }

  cox_fit <- tryCatch(
    survival::coxph(
      survival::Surv(time, status) ~ z,
      data = data,
      ties = "efron",
      x = TRUE,
      y = TRUE,
      model = FALSE
    ),
    warning = function(w) structure(list(warning = conditionMessage(w)), class = "m2_warning"),
    error = function(e) structure(list(error = conditionMessage(e)), class = "m2_error")
  )

  if (inherits(cox_fit, "m2_warning")) {
    warning_message <- cox_fit$warning
    cox_fit2 <- tryCatch(
      suppressWarnings(survival::coxph(
        survival::Surv(time, status) ~ z,
        data = data,
        ties = "efron",
        x = TRUE,
        y = TRUE,
        model = FALSE
      )),
      error = function(e) structure(list(error = e$message), class = "m2_error")
    )
    if (inherits(cox_fit2, "m2_error")) return(m2_failure(cox_fit2$error))
    cox_fit <- cox_fit2
  } else if (inherits(cox_fit, "m2_error")) {
    return(m2_failure(cox_fit$error))
  } else {
    warning_message <- NA_character_
  }

  zph <- tryCatch(
    survival::cox.zph(cox_fit, transform = "identity"),
    error = function(e) structure(list(error = conditionMessage(e)), class = "m2_error")
  )
  if (inherits(zph, "m2_error")) return(m2_failure(zph$error))

  if (is.null(zph$y) || !("z" %in% colnames(zph$y))) {
    return(m2_failure("cox.zph did not return treatment scaled Schoenfeld estimates."))
  }

  smooth_data <- data.frame(
    event_time = as.numeric(zph$x),
    beta = as.numeric(zph$y[, "z"])
  )
  smooth_data <- smooth_data[is.finite(smooth_data$event_time) & is.finite(smooth_data$beta), ]

  unique_events <- length(unique(smooth_data$event_time))
  k <- min(as.integer(basis_dimension), unique_events - 1L)
  if (k < 3L) {
    return(m2_failure("Too few unique event times for the pre-specified smoother."))
  }

  smoother <- tryCatch(
    mgcv::gam(beta ~ s(event_time, k = k, bs = "tp"), data = smooth_data, method = "REML"),
    warning = function(w) structure(list(warning = conditionMessage(w)), class = "m2_warning"),
    error = function(e) structure(list(error = conditionMessage(e)), class = "m2_error")
  )

  if (inherits(smoother, "m2_warning")) {
    smoother_warning <- smoother$warning
    smoother2 <- tryCatch(
      suppressWarnings(mgcv::gam(beta ~ s(event_time, k = k, bs = "tp"), data = smooth_data, method = "REML")),
      error = function(e) structure(list(error = e$message), class = "m2_error")
    )
    if (inherits(smoother2, "m2_error")) return(m2_failure(smoother2$error))
    smoother <- smoother2
  } else if (inherits(smoother, "m2_error")) {
    return(m2_failure(smoother$error))
  } else {
    smoother_warning <- NA_character_
  }

  estimates <- predict_m2_from_smoother(smoother, grid = grid, conf_level = conf_level)

  list(
    method = "M2",
    label = "M2, smoothed scaled Schoenfeld residuals",
    converged = TRUE,
    failure_reason = NA_character_,
    warning = paste(stats::na.omit(c(warning_message, smoother_warning)), collapse = " | "),
    cox_fit = cox_fit,
    zph = zph,
    smoother = smoother,
    smoother_basis_dimension = k,
    event_curve = smooth_data,
    estimates = estimates
  )
}

predict_m2 <- function(fit_object, grid, conf_level = 0.95) {
  if (!isTRUE(fit_object$converged)) return(fit_object)
  predict_m2_from_smoother(fit_object$smoother, grid = grid, conf_level = conf_level)
}

predict_m2_from_smoother <- function(smoother, grid, conf_level = 0.95) {
  stopifnot(all(grid > 0))
  pred <- stats::predict(
    smoother,
    newdata = data.frame(event_time = grid),
    type = "link",
    se.fit = TRUE
  )

  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
  beta <- as.numeric(pred$fit)
  se_beta <- as.numeric(pred$se.fit)
  lower_beta <- beta - zcrit * se_beta
  upper_beta <- beta + zcrit * se_beta

  data.frame(
    method = "M2",
    day = grid,
    log_hr = beta,
    se_log_hr = se_beta,
    ve = 1 - exp(beta),
    lower = 1 - exp(upper_beta),
    upper = 1 - exp(lower_beta),
    conf_level = conf_level
  )
}

m2_failure <- function(reason) {
  list(
    method = "M2",
    label = "M2, smoothed scaled Schoenfeld residuals",
    converged = FALSE,
    failure_reason = reason,
    warning = NA_character_,
    cox_fit = NULL,
    zph = NULL,
    smoother = NULL,
    smoother_basis_dimension = NA_integer_,
    event_curve = NULL,
    estimates = NULL
  )
}

