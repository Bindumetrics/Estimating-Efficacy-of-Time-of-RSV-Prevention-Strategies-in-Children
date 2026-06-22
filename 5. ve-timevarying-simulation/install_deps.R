# install_deps.R
# -----------------------------------------------------------------------------
# Installs every R package this project uses, plus the CmdStan backend that M4
# needs. Safe to re-run: it only installs what is missing (use VE_FORCE=1 to
# reinstall CRAN packages).
#
# Run from project root:  Rscript install_deps.R
# (install_deps.bat finds Rscript on Windows and runs this for you.)
# -----------------------------------------------------------------------------

# ---- CRAN packages used across the project ---------------------------------
cran_pkgs <- c(
  "dplyr",        # evaluate.R / summarise.R aggregation
  "tidyr",        # summarise.R reshaping
  "gt",           # summarise.R HTML tables
  "ggplot2",      # plots.R figures
  "future",       # run_simulation.R parallel backend
  "future.apply", # run_simulation.R parallel map
  "survival",     # M1/M1b/M2 Cox models (recommended pkg; ensure present)
  "mgcv",         # M2 Schoenfeld smoothing (recommended pkg; ensure present)
  "rstpm2",       # M1b supplementary spline Cox
  "posterior"     # M4 posterior summaries
)

repos <- "https://cloud.r-project.org"
force <- tolower(Sys.getenv("VE_FORCE", "")) %in% c("1", "true", "yes")

installed <- rownames(installed.packages())
need <- if (force) cran_pkgs else setdiff(cran_pkgs, installed)

cat("================ R PACKAGE INSTALL ================\n")
if (length(need) == 0) {
  cat("All CRAN packages already installed:\n  ", paste(cran_pkgs, collapse = ", "), "\n", sep = "")
} else {
  cat("Installing CRAN packages:\n  ", paste(need, collapse = ", "), "\n", sep = "")
  install.packages(need, repos = repos)
}

# ---- cmdstanr (not on CRAN; from the Stan R-universe) -----------------------
cat("\n---- cmdstanr (Stan interface for M4) ----\n")
if (!requireNamespace("cmdstanr", quietly = TRUE) || force) {
  install.packages("cmdstanr",
                   repos = c("https://stan-dev.r-universe.dev", repos))
}
have_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
if (!have_cmdstanr) {
  cat("[WARN] cmdstanr did not install. M4 will not run until it does.\n")
}

# ---- C++ toolchain + CmdStan backend ---------------------------------------
if (have_cmdstanr) {
  cat("\n---- CmdStan backend (needs the Rtools C++ toolchain) ----\n")
  toolchain_ok <- tryCatch({
    cmdstanr::check_cmdstan_toolchain(fix = FALSE, quiet = TRUE); TRUE
  }, error = function(e) {
    cat("[WARN] C++ toolchain not ready: ", conditionMessage(e), "\n", sep = "")
    FALSE
  })

  if (!toolchain_ok) {
    cat("\nThe Rtools toolchain is missing or not on PATH.\n",
        "Install the OFFICIAL CRAN Rtools (do NOT use winget here) from:\n",
        "  https://cran.r-project.org/bin/windows/Rtools/\n",
        "Match the Rtools version to your R version, then re-run this installer.\n", sep = "")
  } else {
    cur <- tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE),
                    error = function(e) NA_character_)
    if (is.na(cur)) {
      cat("Toolchain OK. Installing CmdStan (one-time compile, can take minutes)...\n")
      tryCatch(
        cmdstanr::install_cmdstan(cores = max(1L, parallel::detectCores() - 1L),
                                  overwrite = FALSE),
        error = function(e) cat("[WARN] install_cmdstan failed: ",
                                conditionMessage(e), "\n", sep = ""))
    } else {
      cat("CmdStan already installed (version ", cur, ").\n", sep = "")
    }
  }
}

# ---- Final report -----------------------------------------------------------
cat("\n================ VERIFICATION ================\n")
check_one <- function(p) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-14s %s\n", p, if (ok) "OK" else "MISSING"))
  ok
}
all_ok <- all(vapply(c(cran_pkgs, "cmdstanr"), check_one, logical(1)))
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  v <- tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE),
                error = function(e) NA_character_)
  cat(sprintf("  %-14s %s\n", "CmdStan", if (is.na(v)) "NOT INSTALLED" else v))
}
cat("==============================================\n")
if (all_ok) cat("\nAll R packages present. You can now run run_full.bat.\n") else
  cat("\nSome packages are missing - see warnings above.\n")
