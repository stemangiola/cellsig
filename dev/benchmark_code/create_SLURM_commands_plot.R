input_directory = "dev/benchmark_results_multiPC_NH/"
output_directory = "dev/benchmark_results_multiPC_NH/"

tibble(dim = c(2, 4, 10)) %>% 
  mutate(command = sprintf(
    "benchmark_plot_dim%s.pdf:\n\tRscript dev/benchmark_code/produce_plot_from_results.R %s %s %s",
    dim, input_directory, output_directory, dim
    )) %>% 
  
  pull(command) %>% 
  
  purrr::prepend("CATEGORY=yes_no_hierarchy\nMEMORY=60000\nCORES=2\nWALL_TIME=86400") %>% 
  
  write_lines("dev/benchmark_code/benchmark_plot.makeflow")



# if error file needed, add:  > dev/AAA_err.stderr  2>&1 as a suffix to "XXX  > dev/AAA_err.stderr  2>&1" in sprintf