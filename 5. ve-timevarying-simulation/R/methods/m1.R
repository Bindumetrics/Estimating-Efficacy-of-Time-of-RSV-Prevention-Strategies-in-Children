# methods/m1.R
# -----------------------------------------------------------------------------
# M1: Cox model with treatment-by-log-time interaction.
#
# Thesis specification:
#   h(t) = h0(t) * exp(beta1 * Z + beta2 * Z * log(t))
#   VE_h(t) = 1 - exp(beta1 + beta2 * log(t))
#
# This method is intentionally restrictive. It is correctly specified only when
# the true log hazard ratio is linear in log time; otherwise its bias is part of
# the estimand-comparison result, not a coding problem to smooth away.
# -----------------------------------------------------------------------------

fit_m1 <- function(data, grid, conf_level = 0.95) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    return(m1_failure("The survival package is required for M1."))
  }

  required <- c("time", "status", "z")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    return(m1_failure(paste("Missing required columns:", paste(missing, collapse = ", "))))
  }

  if (sum(data$status == 1L, na.rm = TRUE) < 2L) {
    return(m1_failure("Fewer than two events; Cox model is not estimable."))
  }

  fit <- tryCatch(
    survival::coxph(
      survival::Surv(time, status) ~ z + survival::tt(z),
      data = data,
      ties = "efron",
      tt = function(x, t, ...) x * log(pmax(t, .Machine$double.eps)),
      model = FALSE,
      x = FALSE,
      y = FALSE
    ),
    warning = function(w) structure(list(warning = conditionMessage(w)), class = "m1_warning"),
    error = function(e) structure(list(error = conditionMessage(e)), class = "m1_error")
  )

  if (inherits(fit, "m1_warning")) {
    warning_message <- fit$warning
    fit2 <- tryCatch(
      suppressWarnings(survival::coxph(
        survival::Surv(time, status) ~ z + survival::tt(z),
        data = data,
        ties = "efron",
        tt = function(x, t, ...) x * log(pmax(t, .Machine$double.eps)),
        model = FALSE,
        x = FALSE,
        y = FALSE
      )),
      error = function(e) structure(list(error = e$message), class = "m1_error")
    )
    if (inherits(fit2, "m1_error")) return(m1_failure(fit2$error))
    fit <- fit2
  } else if (inherits(fit, "m1_error")) {
    return(m1_failure(fit$error))
  } else {
    warning_message <- NA_character_
  }

  coefs <- stats::coef(fit)
  terms <- identify_m1_terms(names(coefs))
  if (is.null(terms) || any(!is.finite(coefs[terms]))) {
    return(m1_failure("M1 coefficients are missing or non-finite."))
  }

  vc <- tryCatch(stats::vcov(fit)[terms, terms, drop = FALSE], error = function(e) NULL)
  if (is.null(vc) || any(!is.finite(vc))) {
    return(m1_failure("M1 variance-covariance matrix is missing or non-finite."))
  }

  beta <- stats::setNames(coefs[terms], c("z", "log_time_interaction"))
  dimnames(vc) <- list(names(beta), names(beta))

  estimates <- predict_m1_from_coefficients(
    beta = beta,
    vcov = vc,
    grid = grid,
    conf_level = conf_level
  )

  list(
    method = "M1",
    converged = TRUE,
    failure_reason = NA_character_,
    warning = warning_message,
    fit = fit,
    coefficients = beta,
    vcov = vc,
    estimates = estimates
  )
}

identify_m1_terms <- function(coef_names) {
  treatment <- which(coef_names == "z")
  interaction <- grep("tt\\(z\\)|::tt\\(z\\)", coef_names)
  if (length(treatment) != 1L || length(interaction) != 1L) return(NULL)
  c(coef_names[treatment], coef_names[interaction])
}

predict_m1 <- function(fit_object, grid, conf_level = 0.95) {
  if (!isTRUE(fit_object$converged)) return(fit_object)
  predict_m1_from_coefficients(fit_object$coefficients, fit_object$vcov, grid, conf_level)
}

predict_m1_from_coefficients <- function(beta, vcov, grid, conf_level = 0.95) {
  stopifnot(all(grid > 0))
  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
  log_time <- log(grid)
  design <- cbind(1, log_time)

  eta <- as.vector(design %*% beta)
  se_eta <- sqrt(rowSums((design %*% vcov) * design))

  lower_eta <- eta - zcrit * se_eta
  upper_eta <- eta + zcrit * se_eta

  data.frame(
    method = "M1",
    day = grid,
    log_hr = eta,
    se_log_hr = se_eta,
    ve = 1 - exp(eta),
    lower = 1 - exp(upper_eta),
    upper = 1 - exp(lower_eta),
    conf_level = conf_level
  )
}

m1_failure <- function(reason) {
  list(
    method = "M1",
    converged = FALSE,
    failure_reason = reason,
    warning = NA_character_,
    fit = NULL,
    coefficients = NULL,
    vcov = NULL,
    estimates = NULL
  )
}


