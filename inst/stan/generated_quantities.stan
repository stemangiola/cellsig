// functions{
//   	vector lp_reduce_simple( vector global_parameters , vector mus_sigmas , real[] real_data , int[] int_data ) {
// 
//     int G = int_data[1];
// 		int counts[G];
// 		real threshold = -999;
// 		int size_buffer = get_real_buffer_size(mus_sigmas, threshold);
// 		int size_vector = (rows(mus_sigmas)-size_buffer)/2;
// 
// 		if(min(mus_sigmas[1:(size_vector*2)]) == threshold) print("ERROR! The MPI implmentation is buggy");
// 
// 		// Reference / exposure rate
// 		for(g in 1:G)
// 	    counts = neg_binomial_2_log_rng(lambda_log,	1.0 ./ exp( sigma_inv_log )	)
// 	
// 		lp = neg_binomial_2_log_lpmf(
// 			int_data[1:size_vector] |
// 			mus_sigmas[1:size_vector],
// 			1.0 ./ exp( mus_sigmas[size_vector+1:size_vector+size_vector] )
// 		);
// 
// 	 return [lp]';
// 
// 	}
// 	
//   }
data { int<lower=0> G; } 
parameters {

  // Local properties of the data
  vector[G] lambda_log;
  vector[G] sigma_inv_log;

}
generated quantities{

  int counts[G] = neg_binomial_2_log_rng(lambda_log,	1.0 ./ exp( sigma_inv_log )	);
 
}
