# ============================================================
# Bayesian Waning Efficacy Analysis - Medically Attended RSV LRTI
# Fixed version for a SINGLE RDS file: MD_RSV_LRTI.rds
# User directory: D:/Desktop/Tiwonge/Tiwonge
#
# Main fixes in this version:
#   1. Uses the real structure of your RDS object: km_clean, ipd,
#      diagnostics, est, t_max, etc.
#   2. Does NOT assume data_list$placebo or data_list$nirsevimab exist.
#   3. Converts follow-up times to clean, non-negative, whole days.
#   4. Saves all graphs automatically as PNG files, so they appear even
#      when running the script with source() or in a non-interactive session.
#   5. Adds checks that stop early with clear messages if variables are missing.
# ============================================================

# ============================================================
# SECTION 0: USER SETTINGS
# ============================================================

# Your working folder
DATA_DIR <- "D:/Desktop/Tiwonge/Tiwonge"

# Your single R data file for the medically attended RSV LRTI endpoint
DATA_FILE <- file.path(DATA_DIR, "MD_RSV_LRTI.rds")

# Output folder for graphs and saved results
OUT_DIR <- file.path(DATA_DIR, "analysis_outputs_MD_RSV_LRTI_BAYESIAN_ONLY")

# MCMC settings. Increase these for final results.
# For testing, keep them moderate so the script runs quickly.
N_SAMPLES <- 20000
BURN_IN   <- 1000
THIN      <- 20


# Time interval width for aggregated event rates, in days.
# 1 = daily; 7 = weekly. Weekly is often more stable for sparse RSV events.
DT <- 1

# ============================================================
# SECTION 1: SETUP
# ============================================================

setwd(DATA_DIR)

needed_packages <- c("ggplot2", "reshape2", "dplyr", "survival")
for (pkg in needed_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Please install package '", pkg, "' before running this script.")
  }
}

library(ggplot2)
library(reshape2)
library(dplyr)
library(survival)


if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

save_plot <- function(plot_obj, filename, width = 9, height = 6, dpi = 300) {
  print(plot_obj)
  ggplot2::ggsave(
    filename = file.path(OUT_DIR, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
}

message("Working directory: ", DATA_DIR)
message("Output directory: ", OUT_DIR)

# ============================================================
# SECTION 2: LOAD DATA AND BASIC INSPECTION
# ============================================================

data_list <- readRDS(DATA_FILE)

cat("\nTop-level names in RDS object:\n")
print(names(data_list))

cat("\nObject structure summary:\n")
str(data_list, max.level = 2)

# ============================================================
# SECTION 3: HELPER FUNCTIONS FOR ROBUST VARIABLE DETECTION
# ============================================================

clean_days <- function(x) {
  # Converts time to numeric, removes sign errors, and rounds decimals to whole days.
  # This addresses the problem of negative and decimal days.
  x <- suppressWarnings(as.numeric(x))
  x <- abs(x)
  x <- round(x, 0)
  return(x)
}

standardize_arm <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z <- gsub("[_ -]+", "", z)
  out <- ifelse(grepl("plac|control|standard", z), "placebo",
                ifelse(grepl("nir|nirse|treat|active|mab", z), "nirsevimab", z))
  return(out)
}

find_col <- function(df, candidates, required = TRUE, label = "variable") {
  nms <- names(df)
  nms_low <- tolower(nms)
  cand_low <- tolower(candidates)

  hit <- which(nms_low %in% cand_low)
  if (length(hit) > 0) return(nms[hit[1]])

  # Partial matching as fallback
  for (cc in cand_low) {
    hit <- grep(cc, nms_low, fixed = TRUE)
    if (length(hit) > 0) return(nms[hit[1]])
  }

  if (required) {
    stop("Could not find ", label, ". Candidate names tried: ",
         paste(candidates, collapse = ", "),
         "\nAvailable names are: ", paste(nms, collapse = ", "))
  }
  return(NULL)
}

make_survival_from_ipd <- function(ipd) {
  # Builds KM curves from individual participant data if km_clean is missing
  # or if km_clean is not in the expected long format.

  ipd <- as.data.frame(ipd)

  arm_col <- find_col(ipd, c("arm", "group", "treatment", "treat", "trt", "vacc", "vaccine"),
                      label = "treatment arm column in ipd")
  time_col <- find_col(ipd, c("time", "t", "day", "days", "followup", "follow_up", "fu", "ftime"),
                       label = "time column in ipd")
  event_col <- find_col(ipd, c("event", "status", "d", "delta", "case", "infection", "infected", "outcome"),
                        label = "event/status column in ipd")

  tmp <- data.frame(
    arm = standardize_arm(ipd[[arm_col]]),
    t = clean_days(ipd[[time_col]]),
    event = suppressWarnings(as.numeric(ipd[[event_col]]))
  )

  # Convert event coding to 0/1 if needed.
  tmp$event <- ifelse(is.na(tmp$event), 0, tmp$event)
  tmp$event <- ifelse(tmp$event > 0, 1, 0)
  tmp <- tmp[is.finite(tmp$t) & !is.na(tmp$arm), ]

  if (!all(c("placebo", "nirsevimab") %in% unique(tmp$arm))) {
    stop("IPD does not contain both expected arms after standardization. Arms found: ",
         paste(unique(tmp$arm), collapse = ", "))
  }

  fit <- survival::survfit(survival::Surv(t, event) ~ arm, data = tmp)
  sf <- summary(fit)

  km <- data.frame(
    arm = sub("^arm=", "", sf$strata),
    t = clean_days(sf$time),
    Proportion_Free = sf$surv,
    nrisk = sf$n.risk,
    nevent = sf$n.event,
    ncensor = sf$n.censor
  )

  return(km)
}

extract_km_clean <- function(data_list) {
  # Preferred source: data_list$km_clean.
  # Fallback source: create KM from data_list$ipd.

  if (!is.null(data_list$km_clean)) {

    # IMPORTANT FIX:
    # Do NOT start with as.data.frame(data_list$km_clean).
    # In this RDS file, km_clean may be a list whose components have different
    # numbers of rows, for example 67 rows for one arm and 69 rows for another.
    # as.data.frame() tries to combine list elements side-by-side and then fails:
    # "arguments imply differing number of rows: 67, 69".
    # The correct approach is to bind the arms underneath each other using bind_rows().
    if (is.data.frame(data_list$km_clean)) {
      km <- data_list$km_clean
    } else if (is.list(data_list$km_clean)) {
      parts <- lapply(names(data_list$km_clean), function(nm) {
        z <- data_list$km_clean[[nm]]
        if (is.null(z)) return(NULL)
        z <- as.data.frame(z)
        if (nrow(z) == 0) return(NULL)
        if (!"arm" %in% names(z)) z$arm <- nm
        z
      })
      parts <- parts[!vapply(parts, is.null, logical(1))]
      km <- dplyr::bind_rows(parts)
    } else {
      km <- as.data.frame(data_list$km_clean)
    }

    if (nrow(km) > 0) {
      arm_col <- find_col(km, c("arm", "group", "treatment", "treat", "trt"),
                          required = FALSE, label = "arm column in km_clean")
      time_col <- find_col(km, c("t", "time", "day", "days"),
                           required = FALSE, label = "time column in km_clean")
      surv_col <- find_col(km, c("surv", "survival", "Proportion_Free", "proportion_free", "S", "km"),
                           required = FALSE, label = "survival column in km_clean")
      nrisk_col <- find_col(km, c("nrisk", "n.risk", "n_risk", "risk"), required = FALSE)
      nevent_col <- find_col(km, c("nevent", "n.event", "n_event", "events"), required = FALSE)
      ncensor_col <- find_col(km, c("ncensor", "n.censor", "n_censor", "censor"), required = FALSE)

      if (!is.null(arm_col) && !is.null(time_col) && !is.null(surv_col)) {
        out <- data.frame(
          arm = standardize_arm(km[[arm_col]]),
          t = clean_days(km[[time_col]]),
          Proportion_Free = suppressWarnings(as.numeric(km[[surv_col]]))
        )
        if (!is.null(nrisk_col)) out$nrisk <- suppressWarnings(as.numeric(km[[nrisk_col]]))
        if (!is.null(nevent_col)) out$nevent <- suppressWarnings(as.numeric(km[[nevent_col]]))
        if (!is.null(ncensor_col)) out$ncensor <- suppressWarnings(as.numeric(km[[ncensor_col]]))

        out <- out[is.finite(out$t) & is.finite(out$Proportion_Free), ]
        out <- out[out$Proportion_Free >= 0 & out$Proportion_Free <= 1, ]
        if (nrow(out) > 0) return(out)
      }
    }
  }

  if (!is.null(data_list$ipd)) {
    message("km_clean could not be used directly. Reconstructing KM curves from ipd.")
    return(make_survival_from_ipd(data_list$ipd))
  }

  stop("Could not extract KM data. Neither usable km_clean nor ipd was found.")
}

add_zero_anchor <- function(km_df) {
  km_df <- km_df[order(km_df$t, decreasing = FALSE), ]
  km_df <- km_df[!duplicated(km_df$t), ]

  if (nrow(km_df) == 0) stop("KM table has zero rows.")

  # Force clean survival range.
  km_df$Proportion_Free <- pmin(pmax(km_df$Proportion_Free, 0), 1)

  if (!any(km_df$t == 0)) {
    extra <- km_df[1, , drop = FALSE]
    extra$t <- 0
    extra$Proportion_Free <- 1
    if ("nevent" %in% names(extra)) extra$nevent <- 0
    if ("ncensor" %in% names(extra)) extra$ncensor <- 0
    km_df <- rbind(extra, km_df)
  }

  # At day 0 survival must be 1.
  km_df$Proportion_Free[km_df$t == 0] <- 1
  km_df <- km_df[order(km_df$t), ]
  rownames(km_df) <- NULL
  return(km_df)
}

get_n_by_arm <- function(data_list, km_all = NULL) {
  n_by_arm <- NULL

  if (!is.null(data_list$diagnostics) && !is.null(data_list$diagnostics$n_by_arm)) {
    n_by_arm <- data_list$diagnostics$n_by_arm
  }

  if (is.null(n_by_arm) && !is.null(data_list$nrisk)) {
    nr <- as.data.frame(data_list$nrisk)
    arm_col <- find_col(nr, c("arm", "group", "treatment", "treat", "trt"), required = FALSE)
    n_col <- find_col(nr, c("n", "nrisk", "n_risk", "risk"), required = FALSE)
    if (!is.null(arm_col) && !is.null(n_col)) {
      n_by_arm <- tapply(as.numeric(nr[[n_col]]), standardize_arm(nr[[arm_col]]), max, na.rm = TRUE)
    }
  }

  if (is.null(n_by_arm) && !is.null(km_all) && "nrisk" %in% names(km_all)) {
    n_by_arm <- tapply(km_all$nrisk, km_all$arm, max, na.rm = TRUE)
  }

  if (is.null(n_by_arm) && !is.null(data_list$ipd)) {
    ipd <- as.data.frame(data_list$ipd)
    arm_col <- find_col(ipd, c("arm", "group", "treatment", "treat", "trt"), required = FALSE)
    if (!is.null(arm_col)) {
      n_by_arm <- table(standardize_arm(ipd[[arm_col]]))
    }
  }

  if (is.null(n_by_arm)) stop("Could not determine sample size by arm.")

  names(n_by_arm) <- standardize_arm(names(n_by_arm))
  return(n_by_arm)
}

# ============================================================
# SECTION 4: EXTRACT KM DATA CORRECTLY FROM YOUR RDS
# ============================================================

km_all <- extract_km_clean(data_list)

cat("\nKM data extracted. Columns:\n")
print(names(km_all))
cat("\nArms found in KM data:\n")
print(table(km_all$arm, useNA = "ifany"))

# Keep only placebo and nirsevimab arms.
km_all <- km_all[km_all$arm %in% c("placebo", "nirsevimab"), ]

if (!all(c("placebo", "nirsevimab") %in% unique(km_all$arm))) {
  stop("After extracting KM data, both placebo and nirsevimab arms were not found. Arms found: ",
       paste(unique(km_all$arm), collapse = ", "))
}

km_pla <- add_zero_anchor(km_all[km_all$arm == "placebo", ])
km_nir <- add_zero_anchor(km_all[km_all$arm == "nirsevimab", ])

# Keep the core columns used later.
km_pla <- km_pla[, intersect(c("arm", "t", "Proportion_Free", "nrisk", "nevent", "ncensor"), names(km_pla))]
km_nir <- km_nir[, intersect(c("arm", "t", "Proportion_Free", "nrisk", "nevent", "ncensor"), names(km_nir))]

# Follow-up horizon.
T_max <- max(c(km_pla$t, km_nir$t), na.rm = TRUE)
if (!is.null(data_list$t_max)) {
  T_from_object <- clean_days(data_list$t_max)[1]
  if (is.finite(T_from_object) && T_from_object > 0) {
    T_max <- max(T_max, T_from_object, na.rm = TRUE)
  }
}
T_max <- as.integer(ceiling(T_max))

if (!is.finite(T_max) || T_max <= 0) {
  stop("T_max is not positive. Check time variables in km_clean/ipd.")
}

cat("\nFollow-up horizon T_max =", T_max, "days\n")
cat("Placebo KM rows:", nrow(km_pla), "\n")
cat("Nirsevimab KM rows:", nrow(km_nir), "\n")

n_by_arm <- get_n_by_arm(data_list, km_all)
N_pla <- as.numeric(n_by_arm["placebo"])
N_nir <- as.numeric(n_by_arm["nirsevimab"])

cat("\nSample sizes:\n")
cat("  Placebo    =", N_pla, "\n")
cat("  Nirsevimab =", N_nir, "\n")

# ============================================================
# SECTION 5: QUICK KM PLOT
# ============================================================

km_plot_data <- rbind(km_pla, km_nir)

p_km <- ggplot(km_plot_data,
               aes(x = t, y = Proportion_Free, color = arm, group = arm)) +
  geom_step(linewidth = 0.8) +
  coord_cartesian(xlim = c(0, T_max), ylim = c(min(km_plot_data$Proportion_Free, na.rm = TRUE) - 0.005, 1)) +
  labs(title = "Kaplan-Meier Curves - Medically Attended RSV LRTI",
       x = "Days since dosing",
       y = "Proportion free of medically attended RSV LRTI",
       color = "Arm") +
  theme_bw() +
  theme(legend.position = "bottom")

save_plot(p_km, "01_KM_curves.png")

# ============================================================
# SECTION 6: DERIVE EVENT DATA FROM KM CURVES
# ============================================================

km_to_events <- function(km_df, N_total) {
  km_df <- km_df[order(km_df$t), ]
  km_df <- km_df[!duplicated(km_df$t), ]

  events <- data.frame(Type = character(), t = numeric(), n = numeric(),
                       stringsAsFactors = FALSE)

  S_prev <- 1
  at_risk <- N_total

  for (i in seq_len(nrow(km_df))) {
    t_i <- km_df$t[i]
    S_i <- km_df$Proportion_Free[i]

    if (!is.finite(t_i) || !is.finite(S_i)) next
    if (t_i == 0) {
      S_prev <- S_i
      next
    }

    # If actual nevent is available from survfit/summary, use it.
    if ("nevent" %in% names(km_df) && is.finite(km_df$nevent[i]) && km_df$nevent[i] > 0) {
      n_inf <- as.integer(round(km_df$nevent[i]))
    } else {
      delta_S <- S_prev - S_i
      n_inf <- ifelse(delta_S > 1e-10,
                      max(1, round(delta_S / max(S_prev, 1e-12) * at_risk)),
                      0)
    }

    if (n_inf > 0) {
      events <- rbind(events, data.frame(Type = "I", t = t_i, n = n_inf))
      at_risk <- max(0, at_risk - n_inf)
      S_prev <- S_i
    }

    # If censoring count is available, subtract it from at risk.
    if ("ncensor" %in% names(km_df) && is.finite(km_df$ncensor[i]) && km_df$ncensor[i] > 0) {
      n_cen <- as.integer(round(km_df$ncensor[i]))
      events <- rbind(events, data.frame(Type = "C", t = t_i, n = n_cen))
      at_risk <- max(0, at_risk - n_cen)
    }
  }

  # Add remaining people as administratively censored at the end.
  if (at_risk > 0) {
    events <- rbind(events, data.frame(Type = "C", t = max(km_df$t, na.rm = TRUE), n = at_risk))
  }

  events <- events[order(events$t, events$Type), ]
  rownames(events) <- NULL
  return(events)
}

events_pla <- km_to_events(km_pla, N_pla)
events_nir <- km_to_events(km_nir, N_nir)

cat("\nPlacebo events head:\n")
print(head(events_pla, 10))
cat("Total placebo infections:", sum(events_pla$n[events_pla$Type == "I"]), "\n")

cat("\nNirsevimab events head:\n")
print(head(events_nir, 10))
cat("Total nirsevimab infections:", sum(events_nir$n[events_nir$Type == "I"]), "\n")

if (sum(events_pla$n[events_pla$Type == "I"]) == 0 && sum(events_nir$n[events_nir$Type == "I"]) == 0) {
  stop("No infection events were found. Check km_clean or ipd event coding.")
}

# ============================================================
# SECTION 7: BUILD KM ESTIMATOR FROM EVENTS FOR VALIDATION
# ============================================================

KM_est <- function(events) {
  events <- events[order(events$t, events$Type), ]
  events$at_risk_before <- NA_real_

  total_n <- sum(events$n)
  at_risk <- total_n
  S <- 1

  KM_out <- data.frame(t = 0, Proportion_Free = 1)

  for (i in seq_len(nrow(events))) {
    events$at_risk_before[i] <- at_risk

    if (events$Type[i] == "I") {
      KM_out <- rbind(KM_out, data.frame(t = events$t[i], Proportion_Free = S))
      S <- S * (1 - events$n[i] / max(at_risk, 1))
      KM_out <- rbind(KM_out, data.frame(t = events$t[i], Proportion_Free = S))
    } else {
      KM_out <- rbind(KM_out, data.frame(t = events$t[i], Proportion_Free = S))
    }

    at_risk <- max(0, at_risk - events$n[i])
  }

  time_between <- diff(c(0, events$t))
  time_between[time_between < 0] <- 0
  events$pt_elapsed <- cumsum(events$at_risk_before * time_between)

  return(list(Events = events, KM_estimator = KM_out))
}

data_pla_clean <- KM_est(events_pla)
data_nir_clean <- KM_est(events_nir)

recon_plot_data <- rbind(
  data.frame(km_pla[, c("t", "Proportion_Free")], series = "Placebo original"),
  data.frame(km_nir[, c("t", "Proportion_Free")], series = "Nirsevimab original"),
  data.frame(data_pla_clean$KM_estimator, series = "Placebo reconstructed"),
  data.frame(data_nir_clean$KM_estimator, series = "Nirsevimab reconstructed")
)

p_recon <- ggplot(recon_plot_data,
                  aes(x = t, y = Proportion_Free, color = series, group = series)) +
  geom_step(linewidth = 0.7) +
  coord_cartesian(xlim = c(0, T_max), ylim = c(min(recon_plot_data$Proportion_Free, na.rm = TRUE) - 0.005, 1)) +
  labs(title = "KM Curves: Original vs Reconstructed",
       x = "Days since dosing", y = "Proportion free", color = "") +
  theme_bw() +
  theme(legend.position = "bottom")

save_plot(p_recon, "02_KM_original_vs_reconstructed.png")

# ============================================================
# SECTION 8: AGGREGATE INTO TIME INTERVALS
# ============================================================

events_agg <- function(data, ts) {
  events <- data$Events
  output <- data.frame()

  for (i in seq_len(length(ts) - 1)) {
    t_min <- ts[i]
    t_max <- ts[i + 1]

    # Events within interval [t_min, t_max)
    sub <- events[events$t >= t_min & events$t < t_max, ]

    # Approximate at-risk at start of interval.
    past <- events[events$t < t_min, ]
    at_start <- sum(events$n) - ifelse(nrow(past) == 0, 0, sum(past$n))
    at_start <- max(at_start, 0)

    n_inf <- ifelse(nrow(sub) == 0, 0, sum(sub$n[sub$Type == "I"]))
    n_all <- ifelse(nrow(sub) == 0, 0, sum(sub$n))

    # Simple person-time approximation for the interval.
    # Subtract half of events/censors in the interval to account for mid-interval occurrence.
    person_time <- max(0, (at_start - 0.5 * n_all) * (t_max - t_min))

    output <- rbind(output, data.frame(
      t_interval = paste0("[", t_min, ", ", t_max, ")"),
      t = (t_min + t_max) / 2,
      person_time = person_time,
      n = n_inf
    ))
  }

  output$rate <- ifelse(output$person_time > 0, output$n / output$person_time, NA_real_)
  output <- output[is.finite(output$t) & output$person_time > 0, ]
  rownames(output) <- NULL
  return(output)
}

ts <- seq(0, T_max + DT, by = DT)
agg_pla <- events_agg(data_pla_clean, ts)
agg_nir <- events_agg(data_nir_clean, ts)

cat("\nAggregated placebo data head:\n")
print(head(agg_pla))
cat("\nAggregated nirsevimab data head:\n")
print(head(agg_nir))

agg_plot <- rbind(
  data.frame(agg_pla, arm = "Placebo"),
  data.frame(agg_nir, arm = "Nirsevimab")
)

p_rates <- ggplot(agg_plot, aes(x = t, y = rate, color = arm)) +
  geom_point(size = 1.2, alpha = 0.8) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  coord_cartesian(xlim = c(0, T_max), ylim = c(0, NA)) +
  labs(title = "Observed attack rates by interval",
       x = "Days since dosing", y = "Observed attack rate per person-day",
       color = "Arm") +
  theme_bw() +
  theme(legend.position = "bottom")

save_plot(p_rates, "03_observed_attack_rates.png")

# ============================================================
# SECTION 9: BAYESIAN MODEL ONLY
# ============================================================

models <- list(
  "Bayesian" = list(
    name = "Bayesian",
    par_0 = c(2e-4, 7e-4, T_max / 2, 4e-2, -0.95, -1 / (T_max * 2)),
    attack_rate = function(par, t) {
      if (is.null(dim(par))) {
        return(par[1] + par[2] * exp(-(par[4]^2) * (t - par[3])^2 / 2))
      }
      par[,1] + par[,2] * exp(-(par[,4]^2) * (t - par[,3])^2 / 2)
    },
    waning = function(par, t) {
      if (is.null(dim(par))) return(1 - (1 + par[5] * exp(par[6] * t)))
      1 - (1 + par[,5] * exp(par[,6] * t))
    },
    lower = c(0,0,0,0,-1,-1),
    upper = c(1,1,T_max,0.1,0,0),
    log_likelihood = function(par) {
      base_pla <- par[1] + par[2] * exp(-(par[4]^2) * (agg_pla$t - par[3])^2 / 2)
      base_nir <- par[1] + par[2] * exp(-(par[4]^2) * (agg_nir$t - par[3])^2 / 2)

      lambda_pla <- agg_pla$person_time * base_pla
      lambda_nir <- agg_nir$person_time * base_nir *
        (1 + par[5] * exp(par[6] * agg_nir$t))

      sum(
        dpois(agg_pla$n, pmax(lambda_pla,1e-12), log=TRUE),
        dpois(agg_nir$n, pmax(lambda_nir,1e-12), log=TRUE)
      )
    }
  )
)

# ============================================================
# SECTION 10: METROPOLIS MCMC
# ============================================================

log_prior <- function(par, lower, upper) {
  if (any(par < lower) || any(par > upper) || any(!is.finite(par))) return(-Inf)
  sum(dunif(par, lower, upper, log = TRUE))
}

log_posterior <- function(par, model) {
  lp <- log_prior(par, model$lower, model$upper)
  if (!is.finite(lp)) return(-Inf)
  ll <- model$log_likelihood(par)
  if (!is.finite(ll)) return(-Inf)
  lp + ll
}

run_metropolis <- function(model, n_samples = N_SAMPLES, burn_in = BURN_IN, thin = THIN) {
  par <- model$par_0
  step_sd <- pmax(abs(model$par_0) / 5, 1e-6)
  current_lp <- log_posterior(par, model)

  results <- list()
  accept <- 0
  saved <- 0

  for (i in seq_len(n_samples)) {
    par_new <- rnorm(length(par), mean = par, sd = step_sd)
    new_lp <- log_posterior(par_new, model)

    if (is.finite(new_lp) && log(runif(1)) < (new_lp - current_lp)) {
      par <- par_new
      current_lp <- new_lp
      accept <- accept + 1
    }

    if (i > burn_in && i %% thin == 0) {
      saved <- saved + 1
      results[[saved]] <- c(iteration = i, par)
    }
  }

  out <- as.data.frame(do.call(rbind, results))
  names(out) <- c("iteration", paste0("par_", seq_along(model$par_0)))
  attr(out, "acceptance_rate") <- accept / n_samples
  return(out)
}

set.seed(2024)
models_to_run <- which(names(models) == "Bayesian")

for (k in models_to_run) {
  cat("\nRunning Metropolis for model", k, ":", models[[k]]$name, "\n")
  models[[k]]$results <- run_metropolis(models[[k]])
  cat("Acceptance rate:", round(attr(models[[k]]$results, "acceptance_rate"), 3), "\n")
}

# ============================================================
# SECTION 11: PARAMETER SUMMARIES AND TRACE PLOTS
# ============================================================

all_param_summaries <- data.frame()

for (k in models_to_run) {
  results <- models[[k]]$results
  parameter_cols <- grep("^par_", names(results), value = TRUE)

  parameters <- data.frame(
    parameter = parameter_cols,
    mean = sapply(results[parameter_cols], mean, na.rm = TRUE),
    lower_ci = sapply(results[parameter_cols], quantile, probs = 0.025, na.rm = TRUE),
    median = sapply(results[parameter_cols], median, na.rm = TRUE),
    upper_ci = sapply(results[parameter_cols], quantile, probs = 0.975, na.rm = TRUE),
    row.names = NULL
  )
  parameters$model <- models[[k]]$name
  models[[k]]$parameters <- parameters
  all_param_summaries <- rbind(all_param_summaries, parameters)

  results_long <- reshape2::melt(results, id.vars = "iteration", variable.name = "parameter")

  p_trace <- ggplot(results_long, aes(x = iteration, y = value, color = parameter)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~parameter, scales = "free_y") +
    labs(title = paste("Trace plot -", models[[k]]$name),
         x = "Iteration", y = "Parameter value") +
    theme_bw() +
    theme(legend.position = "none")

  save_plot(p_trace, paste0("04_trace_", models[[k]]$name, ".png"), width = 10, height = 7)
}

cat("\nParameter summaries:\n")
print(all_param_summaries)
write.csv(all_param_summaries, file.path(OUT_DIR, "parameter_summaries.csv"), row.names = FALSE)

# ============================================================
# SECTION 12: ATTACK RATE PLOTS
# ============================================================

ts_plot <- sort(unique(agg_pla$t))

for (k in models_to_run) {
  model <- models[[k]]
  results <- model$results
  par_mat <- as.matrix(results[, grep("^par_", names(results)), drop = FALSE])

  attack_rate_pla <- data.frame(t = ts_plot)
  pred_mat_pla <- sapply(ts_plot, function(t_i) model$attack_rate(par_mat, t_i))
  pred_mat_pla <- t(pred_mat_pla)

  attack_rate_pla$mean <- rowMeans(pred_mat_pla, na.rm = TRUE)
  attack_rate_pla$lower_CI <- apply(pred_mat_pla, 1, quantile, 0.025, na.rm = TRUE)
  attack_rate_pla$upper_CI <- apply(pred_mat_pla, 1, quantile, 0.975, na.rm = TRUE)
  models[[k]]$attack_rate_results_pla <- attack_rate_pla

  p_ar_pla <- ggplot() +
    geom_ribbon(data = attack_rate_pla, aes(x = t, ymin = lower_CI, ymax = upper_CI), alpha = 0.2) +
    geom_line(data = attack_rate_pla, aes(x = t, y = mean), linewidth = 0.8) +
    geom_point(data = agg_pla, aes(x = t, y = rate), size = 1, alpha = 0.7) +
    coord_cartesian(xlim = c(0, T_max), ylim = c(0, NA)) +
    labs(x = "Days since dosing", y = "Attack rate per person-day",
         title = paste("Placebo attack rate -", model$name)) +
    theme_bw()

  save_plot(p_ar_pla, paste0("05_placebo_attack_rate_", model$name, ".png"))

  if (!is.null(model$waning)) {
    attack_rate_nir <- data.frame(t = ts_plot)
    pred_mat_nir <- sapply(ts_plot, function(t_i) {
      base <- model$attack_rate(par_mat, t_i)
      eff <- model$waning(par_mat, t_i)
      base * (1 - eff)
    })
    pred_mat_nir <- t(pred_mat_nir)

    attack_rate_nir$mean <- rowMeans(pred_mat_nir, na.rm = TRUE)
    attack_rate_nir$lower_CI <- apply(pred_mat_nir, 1, quantile, 0.025, na.rm = TRUE)
    attack_rate_nir$upper_CI <- apply(pred_mat_nir, 1, quantile, 0.975, na.rm = TRUE)
    models[[k]]$attack_rate_results_nir <- attack_rate_nir

    p_ar_nir <- ggplot() +
      geom_ribbon(data = attack_rate_nir, aes(x = t, ymin = lower_CI, ymax = upper_CI), alpha = 0.2) +
      geom_line(data = attack_rate_nir, aes(x = t, y = mean), linewidth = 0.8) +
      geom_point(data = agg_nir, aes(x = t, y = rate), size = 1, alpha = 0.7) +
      coord_cartesian(xlim = c(0, T_max), ylim = c(0, NA)) +
      labs(x = "Days since dosing", y = "Attack rate per person-day",
           title = paste("Nirsevimab attack rate -", model$name)) +
      theme_bw()

    save_plot(p_ar_nir, paste0("06_nirsevimab_attack_rate_", model$name, ".png"))
  }
}

# ============================================================
# SECTION 13: EFFICACY WANING PLOTS
# ============================================================

treatment_models <- models_to_run[sapply(models, function(x) !is.null(x$waning))]

for (k in treatment_models) {
  model <- models[[k]]
  results <- model$results
  par_mat <- as.matrix(results[, grep("^par_", names(results)), drop = FALSE])

  efficacy <- data.frame(t = ts_plot)
  eff_mat <- sapply(ts_plot, function(t_i) model$waning(par_mat, t_i))
  eff_mat <- t(eff_mat)

  # Keep efficacy within interpretable plotting range.
  efficacy$mean <- pmin(pmax(rowMeans(eff_mat, na.rm = TRUE), 0), 1)
  efficacy$lower_CI <- pmin(pmax(apply(eff_mat, 1, quantile, 0.025, na.rm = TRUE), 0), 1)
  efficacy$upper_CI <- pmin(pmax(apply(eff_mat, 1, quantile, 0.975, na.rm = TRUE), 0), 1)

  models[[k]]$efficacy_results <- efficacy

  p_eff <- ggplot(efficacy, aes(x = t)) +
    geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.2) +
    geom_line(aes(y = mean), linewidth = 0.9) +
    coord_cartesian(xlim = c(0, T_max), ylim = c(0, 1)) +
    labs(x = "Days since dosing",
         y = "Efficacy / protection",
         title = paste("Waning efficacy -", model$name)) +
    theme_bw()

  save_plot(p_eff, paste0("07_waning_efficacy_", model$name, ".png"))
}

# ============================================================
# SECTION 14: BAYESIAN OUTPUT
# ============================================================

bayesian_k <- which(names(models) == "Bayesian")
bayesian_efficacy <- models[[bayesian_k]]$efficacy_results

write.csv(
  bayesian_efficacy,
  file.path(OUT_DIR, "Bayesian_efficacy_over_time.csv"),
  row.names = FALSE
)

# ============================================================
# SECTION 15: SUMMARY TABLE
# ============================================================

get_eff_at_day <- function(eff_df, day) {
  idx <- which.min(abs(eff_df$t - day))
  round(eff_df$mean[idx], 3)
}

summary_table <- data.frame(
  Model = sapply(treatment_models, function(k) models[[k]]$name),
  Mean_Eff_D0 = sapply(treatment_models, function(k) get_eff_at_day(models[[k]]$efficacy_results, 0)),
  Mean_Eff_D30 = sapply(treatment_models, function(k) get_eff_at_day(models[[k]]$efficacy_results, 30)),
  Mean_Eff_D90 = sapply(treatment_models, function(k) get_eff_at_day(models[[k]]$efficacy_results, 90)),
  Mean_Eff_End = sapply(treatment_models, function(k) round(tail(models[[k]]$efficacy_results$mean, 1), 3))
)

cat("\n===== EFFICACY SUMMARY TABLE =====\n")
print(summary_table)
write.csv(summary_table, file.path(OUT_DIR, "efficacy_summary_table.csv"), row.names = FALSE)

# ============================================================
# SECTION 16: SAVE RESULTS
# ============================================================

save(models, agg_pla, agg_nir, data_pla_clean, data_nir_clean,
     km_pla, km_nir, summary_table, T_max, all_param_summaries,
     file = file.path(OUT_DIR, "MD_RSV_LRTI_waning_results_FIXED.RData"))

cat("\nAll results saved in:\n", OUT_DIR, "\n")
cat("Analysis complete.\n")
