// m4_binomial.stan
// -----------------------------------------------------------------------------
// M4: Bayesian binomial interval-count model (thesis Section 4.6, eq. 4.19).
//
// The simulated IPD is binned into a fixed interval grid. For interval k the
// placebo arm contributes a binomial event count with conditional event
// probability pi_p[k] (a free per-interval baseline "discrete hazard"), and the
// vaccine arm contributes a binomial count with
//     pi_v[k] = pi_p[k] * (1 - VE(t_mid[k])).
// VE(t) is one of three candidate functional forms selected by LOOIC, exactly
// the candidate set the thesis uses: exponential, Erlang-3 (Gamma shape 3), and
// power-law. This is NOT a Weibull-baseline survival likelihood; the baseline is
// the saturated per-interval pi_p and the likelihood is binomial throughout.
//
// One model file serves all three forms via `form`. The two waning parameters
// w1, w2 are reused per form; forms that need only one parameter leave w2
// unused in the likelihood (it then samples from its prior and does not affect
// the log-likelihood or LOOIC, which keeps cross-form LOOIC comparable).
// -----------------------------------------------------------------------------

functions {
  // VE(t) for the three candidate forms. Returns a value in (0, ve0] for t >= 0.
  real ve_curve(real t, real ve0, real w1, real w2, int form) {
    if (form == 1) {                 // exponential: ve0 * exp(-rho t)
      return ve0 * exp(-w1 * t);
    } else if (form == 2) {          // Erlang-3 / Gamma shape 3: ve0 * exp(-x)(1 + x + x^2/2)
      real x = w1 * t;
      return ve0 * exp(-x) * (1 + x + 0.5 * square(x));
    } else {                         // power-law (Weibull-shape): ve0 * exp(-(rho t)^gamma)
      return ve0 * exp(-pow(w1 * t, w2));
    }
  }
}

data {
  int<lower=1> K;                    // number of intervals
  array[K] int<lower=0> n_p;         // placebo at-risk at interval start
  array[K] int<lower=0> x_p;         // placebo events in interval
  array[K] int<lower=0> n_v;         // vaccine at-risk at interval start
  array[K] int<lower=0> x_v;         // vaccine events in interval
  vector<lower=0>[K] t_mid;          // interval midpoints (days)
  int<lower=1, upper=3> form;        // 1 = exponential, 2 = erlang3, 3 = power-law
  int<lower=1> G;                    // number of day-grid evaluation points
  vector<lower=0>[G] t_grid;         // day grid for the posterior VE(t) curve
  // Weakly informative prior hyperparameters (set from config in R).
  real<lower=0> ve0_a;               // Beta shape a for ve0 (primary: 2)
  real<lower=0> ve0_b;               // Beta shape b for ve0 (primary: 2)
  real<lower=0> w1_sd;               // half-normal sd for w1 (rho/rate)
  real<lower=0> w2_meanlog;          // lognormal meanlog for w2 (power gamma)
  real<lower=0> w2_sdlog;            // lognormal sdlog for w2
}

parameters {
  real<lower=0, upper=1> ve0;                  // baseline VE at t = 0
  real<lower=0> w1;                            // exp rho / erlang rate / power rho
  real<lower=0> w2;                            // power gamma (unused by forms 1, 2)
  vector<lower=0, upper=1>[K] pi_p;            // per-interval baseline placebo risk
}

model {
  // Priors (weakly informative; NOT centred on the simulated truth).
  ve0 ~ beta(ve0_a, ve0_b);
  w1 ~ normal(0, w1_sd);                       // half-normal via <lower=0>
  w2 ~ lognormal(w2_meanlog, w2_sdlog);        // centred near 1
  pi_p ~ beta(1, 1);                           // flat baseline; binomial dominates

  // Binomial interval-count likelihood.
  for (k in 1:K) {
    real ve_k = ve_curve(t_mid[k], ve0, w1, w2, form);
    real pi_v = pi_p[k] * (1 - ve_k);
    x_p[k] ~ binomial(n_p[k], pi_p[k]);
    x_v[k] ~ binomial(n_v[k], pi_v);
  }
}

generated quantities {
  vector[2 * K] log_lik;             // one entry per binomial observation, for LOO
  vector[G] ve_grid;                 // posterior VE(t) on the day grid

  for (k in 1:K) {
    real ve_k = ve_curve(t_mid[k], ve0, w1, w2, form);
    real pi_v = pi_p[k] * (1 - ve_k);
    log_lik[k]     = binomial_lpmf(x_p[k] | n_p[k], pi_p[k]);
    log_lik[K + k] = binomial_lpmf(x_v[k] | n_v[k], pi_v);
  }
  for (g in 1:G) {
    ve_grid[g] = ve_curve(t_grid[g], ve0, w1, w2, form);
  }
}
