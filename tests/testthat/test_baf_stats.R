
library(covR)
context("baf_stats functions")

ex_file <- system.file("extdata", "sc10_4k.mpileup.gz", package = "covR")

#detach(package:covR, unload=T)

stats <- file_stats(ex_file, verbose=0)
x1 <- read_matrix(ex_file, nrows=stats[2], cols=1:stats[3], verbose=0)


x2 <- baf_stats(x1[,c(1:3, 5, 6)], minq=0)

test_that("baf_stats returns a matrix of appropriate dimensions", {
  expect_equal(class(x2), "matrix")
  expect_equal(ncol(x2), 12)
  expect_equal(nrow(x2), 4000)  
})


x3 <- baf_summary(x2)

test_that("baf_summary returns a vector of appropriate dimension", {
  expect_equal(class(x3), "numeric")
  expect_equal(length(x3), 7)
})


myBed <- matrix(ncol=4, nrow=5)
myBed[,1] <- "Supercontig_1.10"
#myBed[,1] <- "sc12"
myBed[,2] <- c(1, 1001, 1011, 9100, 10600)
myBed[,3] <- c(10, 1020, 1030, 9110, 10607)
myBed[,4] <- paste("gene", 1:5, sep="_")
#myBed


myGenes <- bedify(myBed, x1)

test_that("bedify returns a list of appropriate dimension", {
  expect_equal(class(myGenes), "list")
  expect_equal(length(myGenes), 5)
})


myGenes_stat <- lapply(myGenes, baf_stats)

test_that("baf_stats works with lapply", {
  expect_equal(class(myGenes_stat), "list")
  expect_equal(length(myGenes_stat), 5)  
})


myGenes_sum <- lapply(myGenes_stat, baf_summary)

test_that("baf_sum works with lapply", {
  expect_equal(class(myGenes_sum), "list")
  expect_equal(length(myGenes_sum), 5)
  expect_equal(class(myGenes_sum[[1]]), "numeric")
  expect_equal(length(myGenes_sum[[1]]), 7)
})

