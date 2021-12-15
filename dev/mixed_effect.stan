functions {
    real partial_sum_lpmf(int[] slice_Y,
                        int start, int end,
                        vector mu,
                        vector shape) {
    return neg_binomial_2_log_lupmf(slice_Y | mu[start:end], 1.0 ./ exp(shape[start:end]) );
                               
  }
}
data {
  
  // OBSERVATION LEVEL
  int<lower=1> N;  // total number of observations samples * genes
  int<lower=1> grouping_gene_idx_N[N, 2];  // grouping and gene indicator per observation // first clumn is dataset second is gene
  int Y[N];  // response variable

  // dataset * gene level
  int<lower=1> G; // total genes
  int<lower=1> D; // total datasets, they are unique with genes, so D >> G
  int<lower=1> grouping_gene_idx_D[D, 2];  // grouping and gene indicator per observation // first clumn is dataset second is gene

  // data for group-level effects of ID 1
  //int<lower=1> num_dataset_times_gene;  // number of grouping levels datasets * genes
  //int<lower=1> num_coefficient;  // number of coefficients per level
  
  // Priors
  real assoc_intercept;
  real assoc_slope;
  real assoc_sd_sd;
  real assoc_sd_shape;
  real lambda_mu;
  real lambda_sigma;
  real lambda_skew;
  
  int<lower=1> grainsize;

}
transformed data {
}
parameters {
  vector[G] gene_mean;  // temporary gene_mean for centered predictors
  vector[G] shape;  // shape parameter
  vector<lower=0>[G] gene_sd;  // group-level standard deviations
  vector[D] z_group_level_effect;  // standardized group-level effects

}
transformed parameters {
  vector[D] group_level_effect = gene_sd[grouping_gene_idx_D[,2]] .* z_group_level_effect;
}
model {
  // likelihood including constants

    // initialize linear predictor term
    vector[N] mu = gene_mean[grouping_gene_idx_N[,2]] ;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += group_level_effect[grouping_gene_idx_N[n,1]];
    }
    
  // target += neg_binomial_2_log_lpmf(Y | mu, 1.0 ./ exp(shape[grouping_gene_idx_N[,2]]));

  target += reduce_sum(partial_sum_lupmf, Y, grainsize, mu, shape[grouping_gene_idx_N[,2]]);
                       
  // priors including constants
  //target += student_t_lpdf(gene_mean | 3, 3.3, 3);
  gene_mean ~ skew_normal(lambda_mu, exp(lambda_sigma), lambda_skew);
	
  target += student_t_lpdf(gene_sd | 8, 0, 3) - 1 * student_t_lccdf(0 | 8, 0, 3);
  target += std_normal_lpdf(z_group_level_effect);
  
  // prior association
  shape ~ normal(  gene_mean * assoc_slope + assoc_intercept, assoc_sd_shape);
  //target += normal_lpdf(shape |  0, 3);
  #shape ~  normal(  gene_mean * assoc_slope + assoc_intercept, assoc_sd_sd);
  
}
