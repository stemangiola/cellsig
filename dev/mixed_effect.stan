functions {
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
  #int<lower=1> num_coefficient;  // number of coefficients per level
}
transformed data {
}
parameters {
  vector[G] gene_mean;  // temporary gene_mean for centered predictors
  vector[G] shape;  // shape parameter
  vector[G] gene_sd;  // group-level standard deviations
  vector[D] z_group_level_effect;  // standardized group-level effects
}
transformed parameters {
  vector[D] group_level_effect = exp(gene_sd[grouping_gene_idx_D[,2]]) .* z_group_level_effect;
}
model {
  // likelihood including constants

    // initialize linear predictor term
    vector[N] mu = gene_mean[grouping_gene_idx_N[,2]] ;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += group_level_effect[grouping_gene_idx_N[n,1]];
    }
    
    target += neg_binomial_2_log_lpmf(Y | mu, exp(shape[grouping_gene_idx_N[,2]]));

  // priors including constants
  target += student_t_lpdf(gene_mean | 3, 3.3, 3);
  target += student_t_lpdf(shape | 3, 0, 3);
  target += normal_lpdf(gene_sd | 0, 3)
    - 1 * student_t_lccdf(0 | 3, 0, 3);
  target += std_normal_lpdf(z_group_level_effect);
}
