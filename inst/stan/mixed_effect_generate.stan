functions{

int neg_binomial_2_log_safe_rng(real eta, real phi) {
    real gamma_rate = gamma_rng(phi, phi / exp(eta));
    real log_overflow = 20.7;
    
    if (gamma_rate > exp(log_overflow)) gamma_rate =  exp(log_overflow); // i think this is the max value before overflow but haven't double checked
    if (gamma_rate < exp(-log_overflow)) gamma_rate =  exp(-log_overflow); // i think this is the max value before overflow but haven't double checked

    return poisson_rng(gamma_rate);
  }
  
vector summarise_mean(int[,] grouping_gene_idx_D, vector shape, int G){
  
  int D = num_elements(shape);
  vector[G] shape_mean;
  real placeholder;
  int counter;
  
  for(g in 1:G){
    placeholder = 0.0;
    counter = 0;
    
    for(d in 1:D){
      
      // if gene is g
      if(grouping_gene_idx_D[d,2]==g) {
        placeholder = placeholder + shape[grouping_gene_idx_D[d,1]];
        counter = counter + 1;
      }
      
    }
    
    shape_mean[g] = placeholder/counter;
  }
  
  return(shape_mean);
  
}

}
  data {

  // dataset * gene level
  int<lower=1> G; // total genes
  int<lower=1> D; // total datasets, they are unique with genes, so D >> G
  
  int<lower=1> grouping_gene_idx_D[D, 2];  // grouping and gene indicator per observation // first clumn is dataset second is gene


}
parameters {

  vector[G] gene_mean;  // temporary gene_mean for centered predictors
  vector[G] shape;  // shape parameter
  vector[G] gene_sd;  // group-level standard deviations

}
transformed parameters{
  vector[G] shape_mean = summarise_mean(grouping_gene_idx_D, shape,  G); 

}
generated quantities{

  int Y_gen[G];  // response variable
  real Y_gen_log[G];
  
	for(g in 1:G) {
    	Y_gen[g] = neg_binomial_2_log_safe_rng( normal_rng(gene_mean[g], gene_sd[g]), 1.0 ./ exp(shape_mean[g]));
	}
	
	Y_gen_log = log1p(Y_gen);

}
