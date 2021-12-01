functions {
}
data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  // data for group-level effects of ID 1
  int<lower=1> num_dataset;  // number of grouping levels
  int<lower=1> num_coefficient;  // number of coefficients per level
  int<lower=1> grouping_idx[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] z_group_level_effect_1;
}
transformed data {
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> shape;  // shape parameter
  vector<lower=0>[num_coefficient] sd_1;  // group-level standard deviations
  vector[num_dataset] z_group_level_effect[num_coefficient];  // standardized group-level effects
}
transformed parameters {
  vector[num_dataset] group_level_effect;  // actual group-level effects
  group_level_effect = (sd_1[1] * (z_group_level_effect[1]));
}
model {
  // likelihood including constants

    // initialize linear predictor term
    vector[N] mu = Intercept + rep_vector(0.0, N);
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += group_level_effect[grouping_idx[n]] * z_group_level_effect_1[n];
    }
    target += neg_binomial_2_log_lpmf(Y | mu, shape);

  // priors including constants
  target += student_t_lpdf(Intercept | 3, 3.3, 3);
  target += gamma_lpdf(shape | 0.01, 0.01);
  target += student_t_lpdf(sd_1 | 3, 0, 3)
    - 1 * student_t_lccdf(0 | 3, 0, 3);
  target += std_normal_lpdf(z_group_level_effect[1]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
} 