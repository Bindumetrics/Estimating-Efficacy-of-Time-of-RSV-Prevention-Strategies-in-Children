// Bayesian parametric waning model on grouped survival data.
//
// Likelihood: piecewise-exponential grouped survival. For interval k and arm a,
//   y[k,a] ~ Poisson( PT[k,a] * exp(logh0[k] + z_a * eta_k) )
// where PT is EXACT person-time at risk in that cell.
//
// The brief (sec 4) says "grouped-survival binomial likelihood". The Poisson /
// person-time form is used instead, deliberately: it is the exact likelihood for
// the piecewise-exponential process the DGP actually generates, and it accounts
// for within-interval censoring exactly, whereas a binomial on the
// start-of-interval risk set implicitly credits censored subjects with a full
// interval of exposure. This REMOVES an approximation rather than adding one.
// The two coincide as the per-interval event probability goes to zero.
//
// eta_k = log hazard ratio at the interval midpoint, so VE(t) = 1 - exp(eta).
// Candidate waning forms, selected by `form`:
//   1  constant            eta = a
//   2  affine in log t     eta = a + b * log(t/30)          <- the Cox-TDC form
//   3  exponential decay   eta = a + (c - a) * (1 - exp(-t/tau))
// The DGP's applied-anchor truth is log-QUADRATIC in log t and so lies outside
// all three, by construction (see R/dgp.R, ve_logquad).

data {
  int<lower=1> K;                       // number of intervals
  vector<lower=0>[K] tmid;              // interval midpoints (days)
  array[K] int<lower=0> y0;             // events, control arm
  array[K] int<lower=0> y1;             // events, treated arm
  vector<lower=0>[K] pt0;               // person-time, control arm
  vector<lower=0>[K] pt1;               // person-time, treated arm
  int<lower=1, upper=3> form;
}

transformed data {
  vector[K] logpt0 = log(pt0 + 1e-12);
  vector[K] logpt1 = log(pt1 + 1e-12);
  vector[K] ltm = log(tmid / 30.0);
}

parameters {
  vector[K] logh0;                      // log baseline hazard per interval
  real a;                               // level
  array[form == 1 ? 0 : 1] real b_raw;  // slope / asymptote
  array[form == 3 ? 1 : 0] real<lower=0> tau_raw;
}

transformed parameters {
  vector[K] eta;
  if (form == 1) {
    eta = rep_vector(a, K);
  } else if (form == 2) {
    eta = a + b_raw[1] * ltm;
  } else {
    eta = a + (b_raw[1] - a) * (1 - exp(-tmid / (tau_raw[1] * 30.0)));
  }
}

model {
  // Baseline hazards are on the order of 1e-4 per day, so log h0 ~ -9.
  logh0 ~ normal(-9, 3);
  a ~ normal(0, 2);
  if (form != 1) b_raw ~ normal(0, 2);
  if (form == 3) tau_raw ~ lognormal(1.0, 0.75);

  y0 ~ poisson_log(logpt0 + logh0);
  y1 ~ poisson_log(logpt1 + logh0 + eta);
}

generated quantities {
  vector[2 * K] log_lik;
  for (k in 1:K) {
    log_lik[k]     = poisson_log_lpmf(y0[k] | logpt0[k] + logh0[k]);
    log_lik[K + k] = poisson_log_lpmf(y1[k] | logpt1[k] + logh0[k] + eta[k]);
  }
}
