

context("baf_stats functions")

ex_file <- system.file("extdata", "sc10_4k.mpileup.gz", package = "covR")

library(covR)
stats <- file_stats(ex_file, verbose=0)
x1 <- read_matrix(ex_file, nrows=stats[2], ncols=stats[3], verbose=0)


test_that("baf_stats works",{
  baf_stats(calls=x1[,5], quals=x1[,6], ref=x1[,3], minq=0)
  
})


