
library(covR)
context("Bedify")

ex_file <- system.file("extdata", "sc10_4k.mpileup.gz", package = "covR")
stats <- file_stats(ex_file, verbose=0)
x1 <- read_matrix(ex_file, nrows=stats[2], cols=1:stats[3], verbose=0)


myBed <- matrix(ncol=4, nrow=5)
myBed[,1] <- "sc12"
myBed[,1] <- "Supercontig_1.10"

myBed[,2] <- c(1, 1001, 1011, 9100, 10600)
myBed[,3] <- c(10, 1020, 1030, 9110, 10607)
myBed[,4] <- paste("gene", 1:5, sep="_")


myGenes <- bedify(myBed, x1)
myGenes <- bedify(myBed, x1, fill_missing=1)


test_that("myGenes contains the correct number of genes", {
  expect_equal(length(myGenes), 5)
  
})

