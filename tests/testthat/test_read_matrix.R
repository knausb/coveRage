
library(covR)
context("read_matrix functions")

ex_file <- system.file("extdata", "sc10_4k.mpileup.gz", package = "covR")


stats <- file_stats(ex_file, verbose=0)
#x1 <- read_matrix(ex_file, nrows=stats[2], ncols=stats[3], verbose=0)
#x1 <- read_matrix(ex_file, nrows=stats[2], cols=c(1, 0), verbose=0)
x1 <- read_matrix(ex_file, nrows=stats[2], cols=1:3, verbose=0)
x1 <- read_matrix(ex_file, nrows=stats[2], cols=c(1:3, 53:54), verbose=0)

#x1 <- read_matrix(ex_file, nrows=stats[2], cols=c(1:3, 53:55), verbose=0)

#x1[1:6,]


test_that("Matrix reads in",{
  expect_equal(nrow(x1), 4000)
#  expect_equal(ncol(x1), 54)
})


#x1[,c(4,5)]

