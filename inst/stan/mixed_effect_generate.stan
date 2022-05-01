data {

  // dataset * gene level
  int<lower=1> G; // total genes
  int<lower=1> D; // total datasets, they are unique with genes, so D >> G

}
parameters {

  vector[G] gene_mean;  // temporary gene_mean for centered predictors
  vector[G] shape;  // shape parameter
  vector[G] gene_sd;  // group-level standard deviations

}
generated quantities{

  int Y_gen[G];  // response variable
  real Y_gen_log[G];
  
	for(g in 1:G) {
    	Y_gen[g] = neg_binomial_2_log_rng( normal_rng(gene_mean[g], gene_sd[g]), 1.0 ./ exp(shape[g]));
	}
	
	Y_gen_log = log1p(Y_gen);

}