functions {
  
    vector Q_sum_to_zero_QR(int N) {
    vector [2*N] Q_r;

    for(i in 1:N) {
      Q_r[i] = -sqrt((N-i)/(N-i+1.0));
      Q_r[i+N] = inv_sqrt((N-i) * (N-i+1));
    }
    return Q_r;
  }

  vector sum_to_zero_QR(vector x_raw, vector Q_r) {
    int N = num_elements(x_raw) + 1;
    vector [N] x;
    real x_aux = 0;

    for(i in 1:N-1){
      x[i] = x_aux + x_raw[i] * Q_r[i];
      x_aux = x_aux + x_raw[i] * Q_r[i+N];
    }
    x[N] = x_aux;
    return x;
  }

    real partial_sum_lpmf(int[] slice_Y,
                        int start, int end,
                        vector mu,
                        vector shape, 
                        vector exposure_rate, 
                        int[,] grouping_gene_idx_N) {
                          

    return neg_binomial_2_log_lupmf(
      slice_Y | 
       mu[grouping_gene_idx_N[start:end,1]] +
      
      //exposure
      exposure_rate[grouping_gene_idx_N[start:end,3]],
      1.0 ./ exp(shape[grouping_gene_idx_N[start:end,1]]) 
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
  
  int<lower=1> M; // total multilevel groupings (databases)
  int<lower=1> grouping_gene_idx_D[D, 3];  // grouping and gene indicator per observation // first clumn is dataset second is gene // third column is multilevel grouping database
  vector[S] exposure_rate; // the exposure rate to compensate for sequencing depth

  // data for group-level effects of ID 1
  //int<lower=1> num_dataset_times_gene;  // number of grouping levels datasets * genes
  //int<lower=1> num_coefficient;  // number of coefficients per level
  
  // Hyper Priors
  real assoc_intercept_mean;
  real assoc_slope_mean;
  real<lower=0> assoc_sd_shape_mean;
  real lambda_mu_mean;
  real lambda_sigma_mean;
  real lambda_skew_mean;
  
  // Parallelisation
  int<lower=1> grainsize;
  
  // non centered
  vector[G] gene_mean_offset;

}
transformed data{
 vector[2*M] Q_r = Q_sum_to_zero_QR(M);
  real x_raw_sigma = inv_sqrt(1 - inv(M));
}

parameters {
  vector<offset=gene_mean_offset>[G] gene_mean;  // temporary gene_mean for centered predictors
  
  vector[G] shape_gene;  // shape parameter
  vector[M-1] shape_multilevel_raw;  // shape parameter
  
  vector<lower=0>[G] gene_sd;  // group-level standard deviations
  vector[D] z_group_level_effect;  // standardized group-level effects
  
  real<lower=1> gene_sd_alpha;
  real<lower=1> gene_sd_beta;
   
  //   // Priors
  real<offset = assoc_intercept_mean, multiplier = fabs(assoc_intercept_mean)/5> assoc_intercept;
  real<offset = assoc_slope_mean, multiplier = fabs(assoc_slope_mean)/5> assoc_slope;
  real<lower=0> assoc_sd_shape;
  real<offset = lambda_mu_mean, multiplier = fabs(lambda_mu_mean)/5> lambda_mu;
  real<offset = lambda_sigma_mean, multiplier = fabs(lambda_sigma_mean)/5> lambda_sigma;
  real<offset = lambda_skew_mean, multiplier = fabs(lambda_skew_mean)/5> lambda_skew;

}
transformed parameters {
   vector[M] shape_multilevel =   sum_to_zero_QR(shape_multilevel_raw, Q_r); // shape parameter

  vector[D] group_level_effect = gene_sd[grouping_gene_idx_D[,2]] .* z_group_level_effect;

  vector[D] mu =  
        gene_mean[grouping_gene_idx_D[,2]] +
    
         // add more terms to the linear predictor
        group_level_effect[grouping_gene_idx_D[,1]] ;
        
  vector[D] shape = 
    shape_gene[grouping_gene_idx_D[,2]] +
    shape_multilevel[grouping_gene_idx_D[,3]];
  


}
model {
  // likelihood including constants
  target += reduce_sum(partial_sum_lupmf, Y, grainsize, mu, shape, exposure_rate , grouping_gene_idx_N);

  // priors including constants
  gene_mean ~ skew_normal(lambda_mu, exp(lambda_sigma), lambda_skew);
  gene_sd ~ gamma( gene_sd_alpha, gene_sd_beta); 
  gene_sd_alpha ~ normal(3, 0.5);
  gene_sd_beta ~ normal(3, 0.5);
   
  target += std_normal_lpdf(z_group_level_effect);
  
  // prior association
  shape_gene ~ normal(  gene_mean * assoc_slope + assoc_intercept, assoc_sd_shape);
  shape_multilevel ~ normal(0,1);
  
  // Hyperprior
  assoc_intercept ~ normal(assoc_intercept_mean, fabs(assoc_intercept_mean)/5);
  assoc_slope ~ normal(assoc_slope_mean, fabs(assoc_slope_mean)/5);
  assoc_sd_shape ~ normal(assoc_sd_shape_mean, fabs(assoc_sd_shape_mean)/5);
  lambda_mu ~ normal(lambda_mu_mean, fabs(lambda_mu_mean)/5);
  lambda_sigma ~ normal(lambda_sigma_mean, fabs(lambda_sigma_mean)/5);
  lambda_skew ~ normal(lambda_skew_mean, fabs(lambda_skew_mean)/5);
  
  
}

