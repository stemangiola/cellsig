functions{

}
data {
	// shards
	int<lower=1> shards;

	// Reference matrix inference
	int<lower=0> GM;
	int<lower=0> S;
	int CL; // counts linear size

	// reference counts
 	int<lower=0> counts_linear[CL] ;
	int GM_to_counts_linear[CL] ;
	int S_linear[CL] ;

	// Non-centered param
	real lambda_mu_prior[2];
	real lambda_sigma_prior[2];
	real lambda_skew_prior[2];
	real sigma_intercept_prior[2];

}
parameters {

	// Global properties
	real<offset=lambda_mu_prior[1],multiplier=lambda_mu_prior[2]>lambda_mu;
  real<offset=lambda_sigma_prior[1],multiplier=lambda_sigma_prior[2]> lambda_sigma;
  real<upper=0> lambda_skew;

	// Sigma
	real<upper=0> sigma_slope;
	real<lower=0> sigma_sigma;
  real<offset=sigma_intercept_prior[1],multiplier=sigma_intercept_prior[2]> sigma_intercept;

  // Local properties of the data
  vector[GM] lambda_log;
  vector[GM] sigma_inv_log;
  vector[S-1] exposure_rate_minus_1;


}
transformed parameters{

	vector[S] exposure_rate = append_row(exposure_rate_minus_1, -sum(exposure_rate_minus_1));


}
model {

  // Overall properties of the data
  lambda_mu ~ normal(lambda_mu_prior[1],lambda_mu_prior[2]);
	lambda_sigma ~ normal(lambda_sigma_prior[1],lambda_sigma_prior[2]);
	lambda_skew ~ normal(lambda_skew_prior[1],lambda_skew_prior[2]);

  sigma_intercept ~ normal(0,2);
  sigma_slope ~ normal(0,2);
  sigma_sigma ~ normal(0,2);

	// Exposure
	exposure_rate_minus_1 ~ normal(0,2);

	// Means overdispersion reference
	lambda_log ~ skew_normal(lambda_mu, exp(lambda_sigma), lambda_skew);
	sigma_inv_log ~ normal(sigma_slope * lambda_log + sigma_intercept, sigma_sigma);

  counts_linear ~ neg_binomial_2_log(
    lambda_log[GM_to_counts_linear] + exposure_rate[S_linear],
			1.0 ./ exp( sigma_inv_log[GM_to_counts_linear] )
		);

}
