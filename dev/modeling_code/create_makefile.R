
# library(unixtools)
# dir.create(t <- paste(sprintf("~/.Rtemp/%s", basename(tempdir())), Sys.getpid(), sep='-'), FALSE, TRUE, "0700")
# set.tempdir(t)

library(tidyverse)


local_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/cellsig"
cores = 15

# Create input
sprintf("CATEGORY=create_input\nMEMORY=20024\nCORES=%s", cores) %>%
  
  c(
    dir(sprintf("%s/dev/modeling_results/", local_dir), pattern = "input.rds") %>%
      
      enframe(value = "file") %>%
      mutate(cores = !!cores) %>%
      mutate(command = map2_chr(
        file, cores,
        ~sprintf(
          "dev/modeling_results/%s: dev/modeling_results/%s\n\tRscript dev/modeling_code/run_model.R dev/modeling_results/%s dev/modeling_results/%s %s",
          sprintf("%s_result.rds", basename(.x) %>%  sub("^([^.]*).*", "\\1", .)),
          .x,
          .x,  
          sprintf("%s_result.rds", basename(.x) %>%  sub("^([^.]*).*", "\\1", .)) ,
          .y
        ))
      ) %>%
      pull(command) %>%
      unlist()
  ) %>%
  write_lines("dev/modeling_code/run_model.makeflow") 
