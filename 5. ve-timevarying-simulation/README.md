# Time-varying VE simulation study

This project implements the thesis simulation study in stages. Stage 1 contains the configuration and data-generating mechanism only, so the truth curves, cumulative VE anchors, placebo event rates, randomisation ratios, censoring, and M4 interval grid can be checked before method fitting is added.

## Stage 1 files

- `R/config.R` - all simulation settings and documented provenance notes.
- `R/dgm.R` - calibrated Weibull baseline hazard, hazard-scale VE truths, derived risk-scale VE, inverse-transform event-time simulation, and random loss to follow-up.
- `validate_stage1.R` - creates truth-curve and calibration outputs.
- `tests/test_stage1_dgm.R` - lightweight base-R checks for calibration and simulation behavior.

## Run validation

```r
setwd("C:/Users/bitac/Desktop/ve-timevarying-simulation")
source("validate_stage1.R")
```

or from a terminal:

```bash
Rscript validate_stage1.R
Rscript tests/test_stage1_dgm.R
```

## M4 interval grid

The default interval width is 15 days, matching Chapters 4-5. Change it in `make_sim_config(m4_interval_width = 15L)` when running sensitivity analyses.

## Current scope

This is Stage 1 only. The next stages add M1, M1b, M2, M3, M4, evaluation, orchestration, and reporting.

## Stage 2a: M1

M1 is implemented in `R/methods/m1.R` as the thesis-specified Cox model with treatment-by-log-time interaction:

```r
h(t) = h0(t) * exp(beta1 * Z + beta2 * Z * log(t))
VE_h(t) = 1 - exp(beta1 + beta2 * log(t))
```

Run the single-dataset validation with:

```bash
Rscript validate_m1.R
Rscript tests/test_m1.R
```

The validation output is written to `outputs/m1_validation_*`.

## Stage 2b: M1b supplementary comparator

M1b is implemented in `R/methods/m1b.R` as a supplementary flexible spline Cox comparator on log time. It is **not** one of the four thesis methods and should be labelled as supplementary in all tables and figures.

The current implementation uses `rstpm2::stpm2` with a flexible baseline spline and a time-varying treatment effect. The method remains simulation-only and does not alter the core four-method comparison.

Run validation with:

```bash
Rscript validate_m1b.R
Rscript tests/test_m1b.R
```


## Stage 2c: M2

M2 is implemented in `R/methods/m2.R` as ordinary Cox followed by scaled Schoenfeld residual smoothing. The primary smoother is pre-specified in `config$m2`: an `mgcv` penalised spline fitted by REML, with pointwise intervals from smoother standard errors on the log-HR scale.

Run validation with:

```bash
Rscript validate_m2.R
Rscript tests/test_m2.R
```

## Stage 2d: M3

M3 is implemented in `R/methods/m3.R` directly, without a turnkey package. It fits the Tian/Zucker-Wei kernel-weighted local partial likelihood at each grid day and returns pointwise intervals from the Tian/Zucker-Wei **sandwich** variance `A^{-1} B A^{-1}` (where the meat `B` uses squared kernel weights). Plain inverse information is not used: for a local partial likelihood it is not a consistent variance and produces wildly miscalibrated intervals.

The primary bandwidth and kernel are pre-specified in `config$m3` and should be varied only in sensitivity analyses.

Run validation with:

```bash
Rscript validate_m3.R
Rscript tests/test_m3.R
```

## Stage 2e: M4

M4 is the Bayesian parametric model of Section 4.6, implemented in `R/methods/m4.R` with the Stan model `stan/m4_binomial.stan`. It uses the **binomial interval-count likelihood** (eq. 4.19): the IPD is binned into the fixed interval grid (default 15 days) and, per interval `k`, `pi_v,k = pi_p,k * (1 - VE_r,k)` with `pi_p,k` a free per-interval baseline.

The VE(t) functional form is selected by **LOOIC** from the same three candidates as the thesis — **exponential, Erlang-3 (Gamma shape 3), and power-law**. `fit_m4()` fits all three forms, returns each form's posterior VE(t) curve (for the misspecification sensitivity tier) plus the LOOIC-selected curve (primary).

Primary priors are weakly informative and **not** centred on the truth: `VE0 ~ Beta(2,2)`, waning rate `~ Half-Normal(0, 0.05)`, power shape `~ LogNormal(0, 0.5)`. The matched `Beta(8,2)` prior is sensitivity-tier only.

M4 targets a risk-scale VE; the validation compares it to the DGM-implied **interval-conditional** `VE_r(t)` (not the cumulative VE_r and not the hazard-scale `VE_h`), so a pure scale gap is not mistaken for method bias.

Engine: Stan via `cmdstanr` + CmdStan (needs a C++ toolchain / Rtools on Windows). Run validation with:

```bash
Rscript validate_m4.R
Rscript tests/test_m4.R
```

## Stage 3: evaluation

`R/evaluate.R` harmonises every method onto the common integer-day grid (1..horizon) and scores it against the correct truth **scale**:

- M1, M1b, M2, M3 target hazard-scale `VE_h(t)` and are scored against it.
- M4 targets a risk-scale VE and is scored against the DGM-implied **interval-conditional** `VE_r(t)` (`1 - conditional vaccine risk / conditional placebo risk` per interval). Scoring M4 against `VE_h` would let a pure scale gap look like bias.

Metrics per (method, scenario, n, censoring, grid-day): signed bias, MAE, MSE, RMSE, 95% interval coverage, mean interval width; ISE/IAE over the horizon per replication; and convergence/failure rate as its own outcome (failed estimates are excluded from point summaries but counted in the convergence rate). M3 uses its **pointwise** sandwich intervals for coverage/width. `relative_efficiency()` gives integrated-MSE ratios between methods on the same scale.

**Two comparison layers for RQ4:**

1. **Own-scale (primary, bias-pure):** each method is scored against the estimand it actually targets — M1-M3/M1b vs hazard-scale `VE_h`, M4 vs interval-conditional `VE_r`. This isolates each method's estimation quality.
2. **Common-scale (the four-way head-to-head):** every method is *also* scored against one reference curve (default `VE_h`, "the VE put into the simulation"), giving columns `signed_bias_common`, `mse_common`, etc. `relative_efficiency_common()` ranks all four methods together. M4's risk-scale estimate then carries a small known scale offset; `scale_gap_table(config)` reports `max|VE_h - VE_r|` per scenario so that offset is transparent and not read as estimation error.

**M4 misspecification tier:** `evaluate_methods(..., score_m4_per_form = TRUE)` (the default) expands M4 into the LOOIC-selected curve (`M4`) plus one pseudo-method per candidate form (`M4_exponential`, `M4_erlang3`, `M4_power_law`), so every form is scored against every truth — not just the selected pairing.

`evaluate_methods()` scores one replication; `aggregate_per_day()` / `aggregate_convergence()` / `aggregate_ise()` combine a stack of replications (used by the sweep). Run validation with:

```bash
Rscript validate_evaluate.R    # all methods on one dataset + combined plot
Rscript tests/test_evaluate.R  # deterministic metric-math checks (no Stan)
```

## Stage 4: orchestration

`R/run_simulation.R` runs the tiered sweep. `build_design()` makes the factorial design (all scenarios x n x censoring) at `full_replications`, upgrading headline scenarios at the real-trial-sized anchors to `headline_replications`. `run_simulation()` runs it cell-by-cell (bounding memory), parallel-safe via `future.apply` with **L'Ecuyer-CMRG** streams seeded from `config$meta$rng_seed`; each replication also stores an integer `stan_seed` drawn from its own stream. Failed/non-converged fits yield NA estimates (excluded from point summaries, counted in the convergence rate). It returns and saves aggregated per-day metrics, convergence, ISE, per-replication metadata, and the **M4 LOOIC form-selection table**, as `.rds`.

```bash
Rscript validate_run.R     # tiny SMOKE design incl M4 (reduced sampling) end-to-end
Rscript tests/test_run.R   # design tiering + plumbing + reproducibility (no Stan)
```

For a parallel run: `run_simulation(design, config, parallel = TRUE, workers = N)` — set `controls$m4$parallel_chains = 1` to avoid oversubscription (workers x chains). Use `smoke_design()` / `smoke_controls()` to validate wiring before launching `build_design()` at full replications.

### A note on robust metrics

At very small samples (n = 500, ~15 events) M1's log-time interaction is barely identified, so its `VE(t) = 1 - exp(b1 + b2*log t)` can extrapolate to extreme values at the grid edges. This is a real instability of the naive method (a finding), but it makes the **mean** signed bias and MSE non-robust to a single degenerate replication. `aggregate_per_day()` therefore also reports `median_bias`, `median_abs_error`, and `n_extreme` (replications whose VE left the plausible range). Read mean bias/MSE alongside the robust companions at small n; at n >= 2000 they coincide.

## Stage 5: summaries and figures

`R/summarise.R` turns a sweep result into presentation tables (CSV always; headline tables also as gt HTML): signed bias / RMSE / coverage / mean width by reporting day, ISE, pairwise relative efficiency (within a truth scale), convergence rate, and the M4 LOOIC form-selection table. `R/plots.R` draws per-scenario small-multiples of true VE(t) vs the mean estimated curve with mean bands, each method against the truth on its own scale, with M1b flagged supplementary.

```bash
Rscript validate_summary.R     # tables + figures from the smoke sweep
Rscript tests/test_summarise.R # deterministic table-shape/value checks
```
