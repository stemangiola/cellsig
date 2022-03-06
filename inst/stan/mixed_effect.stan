functions {
    real partial_sum_lpmf(int[] slice_Y,
                        int start, int end,
                        vector gene_mean, 
                        vector group_level_effect,
                        vector shape, 
                        vector exposure_rate, 
                        int[,] grouping_gene_idx_N) {
                          

    return neg_binomial_2_log_lupmf(
      slice_Y | 
       gene_mean[grouping_gene_idx_N[start:end,2]] +
    
      // add more terms to the linear predictor
      group_level_effect[grouping_gene_idx_N[start:end,1]] +
      
      //exposure
      exposure_rate[grouping_gene_idx_N[start:end,3]],
      1.0 ./ exp(shape[grouping_gene_idx_N[start:end,2]]) 
    );
                               
  }
}
data {
  
  // OBSERVATION LEVEL
  int<lower=1> N;  // total number of observations samples * genes
  int<lower=1> grouping_gene_idx_N[N, 3];  // grouping and gene indicator per observation // first clumn is dataset second is gene, third is sample

  int Y[N];  // response variable

  // dataset * gene level
  int<lower=1> G; // total genes
  int<lower=1> D; // total datasets, they are unique with genes, so D >> G
  int<lower=1> S;
  int<lower=1> grouping_gene_idx_D[D, 2];  // grouping and gene indicator per observation // first clumn is dataset second is gene
  vector[S] exposure_rate; // the exposure rate to compensate for sequencing depth

  // data for group-level effects of ID 1
  //int<lower=1> num_dataset_times_gene;  // number of grouping levels datasets * genes
  //int<lower=1> num_coefficient;  // number of coefficients per level
  
  // Priors
  real assoc_intercept;
  real assoc_slope;
  real<lower=0> assoc_sd_shape;
  real lambda_mu;
  real lambda_sigma;
  real lambda_skew;
  
  int<lower=1> grainsize;
  
  // non centered
  vector[G] gene_mean_offset;

}

parameters {
  vector<offset=gene_mean_offset>[G] gene_mean;  // temporary gene_mean for centered predictors
  vector[G] shape;  // shape parameter
  vector<lower=0>[G] gene_sd;  // group-level standard deviations
  vector[D] z_group_level_effect;  // standardized group-level effects
  
  real<lower=1> gene_sd_alpha;
  real<lower=1> gene_sd_beta;
   
  //   // Priors
  // real assoc_intercept;
  // real assoc_slope;
  // real<lower=0> assoc_sd_shape;
  // real lambda_mu;
  // real lambda_sigma;
  // real lambda_skew;

}
transformed parameters {
  vector[D] group_level_effect = gene_sd[grouping_gene_idx_D[,2]] .* z_group_level_effect;
}
model {
  // likelihood including constants

    
  target += reduce_sum(partial_sum_lupmf, Y, grainsize, gene_mean, group_level_effect, shape, exposure_rate , grouping_gene_idx_N);
                       
  // priors including constants
  gene_mean ~ skew_normal(lambda_mu, exp(lambda_sigma), lambda_skew);
	
  gene_sd ~ gamma( gene_sd_alpha, gene_sd_beta); 
  gene_sd_alpha ~ normal(3, 0.5);
  gene_sd_beta ~ normal(3, 0.5);
   
  target += std_normal_lpdf(z_group_level_effect);
  
  // prior association
  shape ~ normal(  gene_mean * assoc_slope + assoc_intercept, assoc_sd_shape);
  //target += normal_lpdf(shape |  0, 3);
  //shape ~  normal(  gene_mean * assoc_slope + assoc_intercept, assoc_sd_sd);
  
}

