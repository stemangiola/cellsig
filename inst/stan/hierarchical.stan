functions{
  vector pcoga_approx(vector x, vector shape, vector rate);
}
data {
  
	// Reference matrix inference
	int<lower=0> G;
	int<lower=0> S;
	int CL; // counts linear size
  int C;

  matrix[S, C+1] prop_template;

	// reference counts
 	int<lower=0> counts_linear[CL] ;

	// Non-centered param
	real lambda_mu_prior[2];
	real lambda_sigma_prior[2];
	real lambda_skew_prior[2];
	real sigma_intercept_prior[2];
	
	// vector[CL] exposure_rate;

}
transformed data{
  vector[3] x = [1,1,1]';
  vector[3] shape = [1,2,3]';
  vector[3] rate = [3,2,1]';
  
  print(pcoga_approx(x, shape, rate));
  
}
parameters {

// 	// Global properties
// 	real<offset=lambda_mu_prior[1],multiplier=lambda_mu_prior[2]>lambda_mu;
//   real<offset=lambda_sigma_prior[1],multiplier=lambda_sigma_prior[2]> lambda_sigma;
//   real<offset=lambda_skew_prior[1],multiplier=lambda_skew_prior[2]> lambda_skew;
// 
// 	// Sigma
// 	real<upper=0> sigma_slope;
// 	real<lower=0> sigma_sigma;
//   real<offset=sigma_intercept_prior[1],multiplier=sigma_intercept_prior[2]> sigma_intercept;

  // Local properties of the data
  matrix[C+1,G] lambda_log;
  matrix[C+1,G] sigma_inv_log;
  
  // vector[C+1] proportion;


}

model {

//   // Overall properties of the data
//   lambda_mu ~ normal(lambda_mu_prior[1],lambda_mu_prior[2]);
// 	lambda_sigma ~ normal(lambda_sigma_prior[1],lambda_sigma_prior[2]);
// 	lambda_skew ~ normal(lambda_skew_prior[1],lambda_skew_prior[2]);

  // sigma_intercept ~ normal(0,2);
  // sigma_slope ~ normal(0,2);
  // sigma_sigma ~ normal(0,2);


	// Means overdispersion reference
	for(c in 1:(C+1)) lambda_log[c] ~ normal(8, 3); // skew_normal(lambda_mu, exp(lambda_sigma), lambda_skew);
	for(c in 1:(C+1)) sigma_inv_log[c] ~ normal(0, 3); //~ normal(sigma_slope * lambda_log + sigma_intercept, sigma_sigma);

  counts_linear ~ neg_binomial_2( to_vector(prop_template * exp(lambda_log)), 1.0 ./ exp(to_vector(sigma_inv_log)));

}

