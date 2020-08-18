functions{
  vector gamma_convoluted(vector shape, vector scale){
    

    
    real mu = sum(shape .* scale);
    real phi =  sum(shape .* square(scale));
    
    real shape_sum = mu^2/phi;
    real scale_sum =  phi/mu;

    return([shape_sum, scale_sum]');
    
  }
}
data {
  
	// Reference matrix inference
	int<lower=0> N;
	int K;

  vector[N] gamma_mix;
  vector<lower=0>[K] shape;

}
parameters {


  vector<lower=0>[K] rate;


}

model {

  vector[2] shape_scale = gamma_convoluted(shape, 1.0 ./ (rate/100));
  shape ~ normal(0,2.5);
  rate ~ normal(0,2.5);
  
  // print(shape, rate);
  // print(shape_rate);
  gamma_mix ~ gamma(shape_scale[1], 1.0 ./shape_scale[2]);

}

