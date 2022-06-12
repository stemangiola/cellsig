source("dev/jian_R_files/test_main_function.R")
library(testthat)
test_that("Main function DE approach", {
  expect_output(str(sig_de_rank), "$ signature", fixed = TRUE)
})
