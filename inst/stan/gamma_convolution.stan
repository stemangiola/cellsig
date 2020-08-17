functions{
  vector gamma_2_convoluted(vector mu, vector phi){
    

    return([sum(mu), sum(phi)]');
    
  }
  
  real gamma_2_lpdf(vector y, real mu, real phi){
    
     real alpha = mu .* mu / phi; 
     real beta = mu / phi;
     
     return( gamma_lpdf(y | alpha, beta) );
  
  }
}
data {
  
	// Reference matrix inference
	int<lower=0> N;
	int K;

  vector[N] gamma_mix;

}
parameters {

  vector[K] log_mu;
  vector<lower=0>[K] phi;


}

model {

  vector[2] mu_phi = gamma_2_convoluted(exp(log_mu), exp(log_mu) + (square(exp(log_mu)) ./ phi));
  log_mu ~ student_t(3,0,2.5);
  phi ~ student_t(3,0,2.5);
  
  print(mu_phi);
  // print(shape_rate);
  //gamma_mix ~ gamma_2(exp(log_mu[1]), exp(log_mu[1]) + (exp(log_mu[1])^2/phi[1]) );
  gamma_mix ~ gamma_2(mu_phi[1], mu_phi[2]);

}

