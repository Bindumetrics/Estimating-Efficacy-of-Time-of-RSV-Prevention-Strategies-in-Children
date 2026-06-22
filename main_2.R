# ============================================================
# RUN ALL METHODS ON MULTIPLE IPD DATASETS (CLEAN PIPELINE)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(survival)
})

# ---- source toolkits
source("methods/dtest.R")      # Durham toolkit (durham_* + plot_durham_ve + bootstrap helpers)
source("methods/tdc.R")        # TDC Cox toolkit
source("methods/logistic.R")   # discrete-time logistic toolkit
source("methods/bayesian.R")   # Stan toolkit (fit_bayes_stan, pp_make, extract_ve_curve_stan, plot_ve_stan)

set.seed(42)

# ============================================================
# 1) DATASETS
# ============================================================
datasets <- tibble::tribble(
  ~study,     ~endpoint,             ~file,
  "Simoes",   "RSV_Hospitalization",  "RSV_Hospitalization.rds",
  "Simoes",   "MD_RSV_LRTI",          "MD_RSV_LRTI.rds",
  "Kampmann", "Severe_RSV_LRTI",      "Severe_RSV_LRTI.rds"
)

data_dir <- "data/reconstructed_ipd"
out_dir  <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 2) HELPERS
# ============================================================

qc_ipd <- function(ipd) {
  stopifnot(all(c("time", "status", "treat") %in% names(ipd)))
  tibble(
    n = nrow(ipd),
    events = sum(ipd$status == 1, na.rm = TRUE),
    time_min = min(ipd$time, na.rm = TRUE),
    time_max = max(ipd$time, na.rm = TRUE),
    n_missing = sum(is.na(ipd$time) | is.na(ipd$status) | is.na(ipd$treat))
  )
}

detect_arms <- function(ipd, ref_hint = "placebo") {
  tr <- unique(as.character(ipd$treat))
  tr <- sort(tr)

  # prefer placebo as reference if present
  if (ref_hint %in% tr) {
    ref <- ref_hint
    vax <- setdiff(tr, ref)
    if (length(vax) != 1) stop("Could not uniquely identify vaccine arm. Levels: ", paste(tr, collapse=", "))
    vax <- vax[1]
  } else {
    # otherwise pick first as ref and second as vax
    if (length(tr) != 2) stop("treat must have exactly 2 levels. Levels: ", paste(tr, collapse=", "))
    ref <- tr[1]
    vax <- tr[2]
  }
  list(ref_level = ref, vax_level = vax)
}

safe_run <- function(expr) {
  tryCatch(expr, error = function(e) e)
}

# extract a simple summary row from each method output
summarise_method <- function(method, res, meta) {
  # meta includes study, endpoint, file, arms, qc
  if (inherits(res, "error")) {
    return(tibble(
      study = meta$study,
      endpoint = meta$endpoint,
      file = meta$file,
      method = method,
      ref_level = meta$ref_level,
      vax_level = meta$vax_level,
      n = meta$qc$n,
      events = meta$qc$events,
      status = "ERROR",
      message = conditionMessage(res)
    ))
  }

  # Method-specific summaries (keep it minimal + consistent)
  if (method == "Durham") {
    return(tibble(
      study = meta$study, endpoint = meta$endpoint, file = meta$file,
      method = method,
      ref_level = meta$ref_level, vax_level = meta$vax_level,
      n = meta$qc$n, events = meta$qc$events,
      status = "OK",
      hr_overall = res$cox$hr,
      ve_overall = res$cox$ve,
      ph_p = res$ph$p_treat,
      aic = AIC(res$cox$fit)
    ))
  }

  if (method == "TDC_Cox") {
    fit_obj <- res$fit_obj
    b0 <- as.numeric(fit_obj$coef["b0"])
    b1 <- as.numeric(fit_obj$coef["b1"])
    ve150 <- tdc_predict_logt(fit_obj, t_grid = 150)$VE_t
    return(tibble(
      study = meta$study, endpoint = meta$endpoint, file = meta$file,
      method = method,
      ref_level = meta$ref_level, vax_level = meta$vax_level,
      n = meta$qc$n, events = meta$qc$events,
      status = "OK",
      b0 = b0, b1 = b1,
      ve_150 = as.numeric(ve150),
      aic = AIC(fit_obj$fit)
    ))
  }

  if (method == "Logit_TDC") {
    fit <- res$fit_obj$fit
    co <- coef(fit)
    # may or may not exist depending on baseline option
    b0 <- unname(co["vax"])
    b1 <- unname(co[grep("vax:logt", names(co))][1])
    return(tibble(
      study = meta$study, endpoint = meta$endpoint, file = meta$file,
      method = method,
      ref_level = meta$ref_level, vax_level = meta$vax_level,
      n = meta$qc$n, events = meta$qc$events,
      status = "OK",
      b0 = as.numeric(b0),
      b1 = as.numeric(b1),
      aic = AIC(fit)
    ))
  }

  if (method == "Bayes_Stan") {
    # you can extend this with Rhat/ESS summaries later
    return(tibble(
      study = meta$study, endpoint = meta$endpoint, file = meta$file,
      method = method,
      ref_level = meta$ref_level, vax_level = meta$vax_level,
      n = meta$qc$n, events = meta$qc$events,
      status = "OK"
    ))
  }

  tibble(
    study = meta$study, endpoint = meta$endpoint, file = meta$file,
    method = method,
    status = "OK"
  )
}

# ============================================================
# 3) RUNNER FOR ONE DATASET
# ============================================================
run_one_dataset <- function(study, endpoint, file,
                            B = 500, seed = 42,
                            eps = 1e-3,
                            dt = 7,          # weekly for logistic + stan preprocessing
                            durham_span = 0.75,
                            bayes_iter = 2000, bayes_warmup = 1000, bayes_chains = 4) {

  path <- file.path(data_dir, file)
  bundle <- readRDS(path)
  ipd <- bundle$ipd

  qc <- qc_ipd(ipd)
  arms <- detect_arms(ipd, ref_hint = "placebo")

  meta <- list(
    study = study, endpoint = endpoint, file = file,
    qc = qc,
    ref_level = arms$ref_level,
    vax_level = arms$vax_level
  )

  # ---- Durham
  durham_res <- safe_run({
    durham_waning_analysis(
      data = ipd,
      time_col = "time", status_col = "status", treat_col = "treat",
      ref_level = arms$ref_level, vax_level = arms$vax_level,
      smooth = "loess", span = durham_span,
      zph_transform = "identity",
      B = B, seed = seed
    )
  })

  # ---- TDC Cox
  tdc_res <- safe_run({
    tdc_logt_analysis(
      data = ipd,
      time_col = "time", status_col = "status", treat_col = "treat",
      ref_level = arms$ref_level, vax_level = arms$vax_level,
      B = B, seed = seed, eps = eps
    )
  })

  # ---- Logistic discrete-time TDC (weekly recommended)
  logit_res <- safe_run({
    logit_tdc_analysis(
      data = ipd,
      time_col="time", status_col="status", treat_col="treat",
      ref_level = arms$ref_level, vax_level = arms$vax_level,
      dt = dt,
      baseline = "factor",
      B = min(300, B), seed = seed
    )
  })

  # ---- Bayesian Stan (weekly)
  bayes_res <- safe_run({
    # requires rstan configured
    pp_obj <- pp_make(ipd, dt = dt, ref_level = arms$ref_level, vax_level = arms$vax_level)
    fit <- fit_bayes_stan(
      pp_obj,
      df_baseline = 8,
      df_waning = 6,
      iter = bayes_iter,
      warmup = bayes_warmup,
      chains = bayes_chains,
      seed = seed,
      refresh = 100
    )
    curve <- extract_ve_curve_stan(fit)
    list(fit = fit, curve = curve)
  })

  # ---- Summaries
  summary_tbl <- bind_rows(
    summarise_method("Durham", durham_res, meta),
    summarise_method("TDC_Cox", tdc_res, meta),
    summarise_method("Logit_TDC", logit_res, meta),
    summarise_method("Bayes_Stan", bayes_res, meta)
  )

  # ---- Curves (standardize column names)
  curves <- list(
    Durham = if (!inherits(durham_res, "error")) durham_res$curve %>% select(time, VE_t, VE_lower, VE_upper) else NULL,
    TDC_Cox = if (!inherits(tdc_res, "error")) tdc_res$curve %>% select(time, VE_t, VE_lower, VE_upper) else NULL,
    Logit_TDC = if (!inherits(logit_res, "error")) logit_res$curve %>% select(time, VE_t, VE_lower, VE_upper) else NULL,
    Bayes_Stan = if (!inherits(bayes_res, "error")) bayes_res$curve %>% rename(VE_t = VE_median, VE_lower = VE_lower, VE_upper = VE_upper) %>% select(time, VE_t, VE_lower, VE_upper) else NULL
  )

  list(
    meta = meta,
    qc = qc,
    summary = summary_tbl,
    results = list(
      durham = durham_res,
      tdc = tdc_res,
      logit = logit_res,
      bayes = bayes_res
    ),
    curves = curves
  )
}

# ============================================================
# 4) RUN ALL DATASETS
# ============================================================

all_out <- vector("list", nrow(datasets))

for (i in seq_len(nrow(datasets))) {
  ds <- datasets[i, ]

  message("\n============================================================")
  message("Running: ", ds$study, " | ", ds$endpoint, " | ", ds$file)
  message("============================================================")

  out <- run_one_dataset(
    study = ds$study,
    endpoint = ds$endpoint,
    file = ds$file,
    B = 500, seed = 42,
    eps = 1e-3,
    dt = 7,
    bayes_iter = 2000, bayes_warmup = 1000, bayes_chains = 4
  )

  all_out[[i]] <- out

  # Save per-dataset objects
  tag <- paste(ds$study, ds$endpoint, sep = "_")
  saveRDS(out, file.path(out_dir, paste0(tag, "_all_methods.rds")))

  # Save summaries + curves
  readr::write_csv(out$summary, file.path(out_dir, paste0(tag, "_summary.csv")))

  for (nm in names(out$curves)) {
    if (!is.null(out$curves[[nm]])) {
      readr::write_csv(out$curves[[nm]], file.path(out_dir, paste0(tag, "_curve_", nm, ".csv")))
    }
  }
}

# Combine summaries across datasets
summary_all <- bind_rows(lapply(all_out, `[[`, "summary"))
readr::write_csv(summary_all, file.path(out_dir, "ALL_summary.csv"))

summary_all




# Example: load your reconstructed IPD bundles (adapt paths to your project)
ipd_severe <- readRDS("data/reconstructed_ipd/RSV_Hospitalization.rds")$ipd
ipd_any    <- readRDS("data/reconstructed_ipd/Any_RSV_LRTI.rds")$ipd   # or whatever you named it

# If your arm labels differ, set them to match your data
arm_levels <- c("placebo", "nirsevimab")

fig2a <- plot_km_cuminc(
  ipd_severe,
  time_var   = "time",
  status_var = "status",
  arm_var    = "treat",
  arm_levels = arm_levels,
  max_day    = 180,
  title      = "Medically Attended Severe RSV-Associated LRTI (Cumulative Incidence)"
)

fig2b <- plot_km_cuminc(
  ipd_any,
  time_var   = "time",
  status_var = "status",
  arm_var    = "treat",
  arm_levels = arm_levels,
  max_day    = 180,
  title      = "Medically Attended RSV-Associated LRTI (Any) (Cumulative Incidence)"
)

fig2a$plot
fig2b$plot

# Risk tables (0,30,60,90,120,150,180) like the paper
fig2a$risk
fig2b$risk


library(survival)
library(ggplot2)
library(dplyr)
library(survminer)

ipd <- ipd_severe

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
  annotate("text", x = 95, y = 0.905, label = hr_txt, hjust = 0, size = 3.5)

simoes_plot()