# methods/m3.R
# -----------------------------------------------------------------------------
# M3: Tian / Zucker-Wei kernel-weighted local partial likelihood.
#
# At each grid day t0, beta(t0) is estimated by solving a kernel-weighted Cox
# partial likelihood score equation. Event contributions are weighted by their
# distance from t0, while risk sets remain the usual Cox risk sets.
#
# For a single treatment indicator Z:
#   U(beta; t0) = sum_i d_i K_h(T_i - t0) [Z_i - E_beta{Z | R(T_i)}]
#
# Pointwise variance is the Tian/Zucker-Wei sandwich, NOT plain inverse
# information. For a local (kernel-weighted) partial likelihood the model-based
# inverse information is not a consistent variance: it scales with the bandwidth
# and is badly miscalibrated. The correct pointwise variance is
#   Var(beta_hat(t0)) = A(t0)^{-1} B(t0) A(t0)^{-1},
#   A(t0) = sum_i d_i K_h(T_i - t0)   V_i,    (weighted information)
#   B(t0) = sum_i d_i K_h(T_i - t0)^2 V_i,    (meat: SQUARED kernel weights)
# with V_i the risk-set variance of Z at T_i. For a scalar beta this is
# se = sqrt(B) / A. The sandwich is invariant to the 1/h kernel normalisation,
# so the point estimate (score / A) is unchanged by this correction.
#
# A simultaneous Tian band can be added later descriptively, but pointwise
# intervals are the like-for-like metric required for coverage comparisons.
# -----------------------------------------------------------------------------

fit_m3 <- function(data, grid, bandwidth_days = 30, kernel = "gaussian", conf_level = 0.95) {
  required <- c("time", "status", "z")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    return(m3_failure(paste("Missing required columns:", paste(missing, collapse = ", "))))
  }
  if (bandwidth_days <= 0) {
    return(m3_failure("bandwidth_days must be positive."))
  }

  data <- data[, required]
  data <- data[is.finite(data$time) & is.finite(data$status) & is.finite(data$z), ]
  data$z <- as.numeric(data$z)
  data$status <- as.integer(data$status)

  if (sum(data$status == 1L) < 6L) {
    return(m3_failure("Fewer than six events; local partial likelihood is not stable."))
  }
  if (length(unique(data$z)) < 2L) {
    return(m3_failure("Both treatment arms are required for M3."))
  }

  event_data <- data[data$status == 1L, , drop = FALSE]
  event_times <- event_data$time
  event_z <- event_data$z

  estimates <- do.call(rbind, lapply(grid, function(t0) {
    fit_m3_at_time(
      t0 = t0,
      data = data,
      event_times = event_times,
      event_z = event_z,
      bandwidth_days = bandwidth_days,
      kernel = kernel,
      conf_level = conf_level
    )
  }))

  converged <- all(estimates$local_converged)
  failure_reason <- if (converged) NA_character_ else "One or more grid-day local likelihood fits failed."

  list(
    method = "M3",
    label = "M3, Tian/Zucker-Wei kernel-weighted local partial likelihood",
    converged = converged,
    failure_reason = failure_reason,
    warning = NA_character_,
    bandwidth_days = bandwidth_days,
    kernel = kernel,
    estimates = estimates
  )
}

fit_m3_at_time <- function(t0, data, event_times, event_z, bandwidth_days, kernel, conf_level) {
  weights <- kernel_weights((event_times - t0) / bandwidth_days, kernel = kernel) / bandwidth_days
  active <- is.finite(weights) & weights > max(weights, na.rm = TRUE) * 1e-8

  if (sum(active) < 3L || length(unique(event_z[active])) < 1L) {
    return(m3_day_failure(t0, conf_level, "Too few locally weighted events."))
  }

  event_times_active <- event_times[active]
  event_z_active <- event_z[active]
  weights_active <- weights[active]

  objective <- function(beta) m3_score_info(beta, data, event_times_active, event_z_active, weights_active)

  beta <- 0
  converged <- FALSE
  for (iter in seq_len(50L)) {
    si <- objective(beta)
    if (!is.finite(si$score) || !is.finite(si$info) || si$info <= 0) break
    step <- si$score / si$info
    step <- max(min(step, 2), -2)
    beta_new <- beta + step
    if (!is.finite(beta_new)) break
    if (abs(beta_new - beta) < 1e-8) {
      beta <- beta_new
      converged <- TRUE
      break
    }
    beta <- beta_new
  }

  si <- objective(beta)
  if (!converged && is.finite(si$score) && is.finite(si$info) && si$info > 0 && abs(si$score) < 1e-5) {
    converged <- TRUE
  }
  if (!converged || !is.finite(si$info) || si$info <= 0) {
    return(m3_day_failure(t0, conf_level, "Local Newton solve failed."))
  }

  # Tian/Zucker-Wei sandwich pointwise standard error: A^{-1} B A^{-1}.
  if (!is.finite(si$meat) || si$meat <= 0) {
    return(m3_day_failure(t0, conf_level, "Non-finite sandwich variance."))
  }
  se <- sqrt(si$meat) / si$info
  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
  lower_beta <- beta - zcrit * se
  upper_beta <- beta + zcrit * se

  data.frame(
    method = "M3",
    day = t0,
    log_hr = beta,
    se_log_hr = se,
    ve = 1 - exp(beta),
    lower = 1 - exp(upper_beta),
    upper = 1 - exp(lower_beta),
    conf_level = conf_level,
    local_events = sum(active),
    effective_events = sum(weights_active) ^ 2 / sum(weights_active ^ 2),
    local_converged = TRUE,
    local_failure_reason = NA_character_
  )
}

m3_score_info <- function(beta, data, event_times, event_z, weights) {
  score <- 0
  info <- 0   # A(t0) = sum w_j V_j        (weighted information; the bread)
  meat <- 0   # B(t0) = sum w_j^2 V_j      (sandwich meat; squared weights)
  exp_beta_z <- exp(beta * data$z)

  for (j in seq_along(event_times)) {
    risk <- data$time >= event_times[j]
    s0 <- sum(exp_beta_z[risk])
    s1 <- sum(data$z[risk] * exp_beta_z[risk])
    s2 <- sum((data$z[risk] ^ 2) * exp_beta_z[risk])

    if (s0 <= 0) next
    mean_z <- s1 / s0
    var_z <- s2 / s0 - mean_z ^ 2

    score <- score + weights[j] * (event_z[j] - mean_z)
    info <- info + weights[j] * var_z
    meat <- meat + (weights[j] ^ 2) * var_z
  }

  list(score = score, info = info, meat = meat)
}

kernel_weights <- function(u, kernel = "gaussian") {
  switch(
    kernel,
    gaussian = stats::dnorm(u),
    epanechnikov = ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0),
    stop("Unknown M3 kernel: ", kernel)
  )
}

m3_day_failure <- function(t0, conf_level, reason) {
  data.frame(
    method = "M3",
    day = t0,
    log_hr = NA_real_,
    se_log_hr = NA_real_,
    ve = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    conf_level = conf_level,
    local_events = NA_integer_,
    effective_events = NA_real_,
    local_converged = FALSE,
    local_failure_reason = reason
  )
}

predict_m3 <- function(fit_object, grid, conf_level = 0.95) {
  if (is.null(fit_object$estimates)) return(fit_object)
  out <- fit_object$estimates[match(grid, fit_object$estimates$day), ]
  out$conf_level <- conf_level
  out
}

m3_failure <- function(reason) {
  list(
    method = "M3",
    label = "M3, Tian/Zucker-Wei kernel-weighted local partial likelihood",
    converged = FALSE,
    failure_reason = reason,
    warning = NA_character_,
    bandwidth_days = NA_real_,
    kernel = NA_character_,
    estimates = NULL
  )
}
