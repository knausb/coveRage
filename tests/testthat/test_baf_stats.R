
library(coveRage)
#library(testthat)
context("baf_stats functions")

ex_file <- system.file("extdata", "sc10_4k.mpileup.gz", package = "coveRage")

#detach(package:covR, unload=T)

stats <- file_stats(ex_file, verbose=0)
x1 <- read_matrix(ex_file, nrows=stats[2], cols=1:stats[3], verbose=0)


x2 <- baf_stats(x1[,c(1:3, 5, 6)], minq=0)

test_that("baf_stats returns a matrix of appropriate dimensions", {
  expect_equal(class(x2), "matrix")
  expect_equal(ncol(x2), 12)
  expect_equal(nrow(x2), 4000)  
})


test_that("baf_stats minq filters calls", {
  x2 <- baf_stats(x1[,c(1:3, 5, 6)], minq=0)
  expect_match(as.character(sum(rowSums(x2[,-1]))), "^1833[[:digit:]]{2}")
  x2 <- baf_stats(x1[,c(1:3, 5, 6)], minq=20)
  expect_match(as.character(sum(rowSums(x2[,-1]))), "^1823[[:digit:]]{2}")
  x2 <- baf_stats(x1[,c(1:3, 5, 6)], minq=30)
  expect_match(as.character(sum(rowSums(x2[,-1]))), "^1765[[:digit:]]{2}")
  x2 <- baf_stats(x1[,c(1:3, 5, 6)], minq=40)
  expect_match(as.character(sum(rowSums(x2[,-1]))), "^841[[:digit:]]{2}")
  
# sum(rowSums(x2[,-1]))
# baf_plot(x2)
})

##### ##### ##### ##### #####
# Test that counts are accurate.


ex_file <- system.file("extdata", "example_variants.mpileup.gz", package = "coveRage")
stats <- file_stats(ex_file, verbose=0)
pile <- read_matrix(ex_file, nrows=stats[2], cols=1:stats[3], verbose=0)

baf <- baf_stats(pile)

test_that("baf_stats counts correctly", {
  # A
  expect_equal(as.numeric(baf[1,2]), 5)
  expect_equal(sum(as.numeric(baf[1,-c(1,2)])), 0)
  # C
  expect_equal(as.numeric(baf[2,3]), 5)
  expect_equal(sum(as.numeric(baf[2,-c(1,3)])), 0)
  # G
  expect_equal(as.numeric(baf[3,4]), 6)
  expect_equal(sum(as.numeric(baf[3,-c(1,4)])), 0)
  # T
  expect_equal(as.numeric(baf[4,5]), 7)
  expect_equal(sum(as.numeric(baf[4,-c(1,5)])), 0)

  # Aa
  expect_equal(as.numeric(baf[5,2]), 9)
  expect_equal(as.numeric(baf[5,8]), 2)
  expect_equal(sum(as.numeric(baf[5,-c(1,2,8)])), 0)
  # Cc
  expect_equal(as.numeric(baf[6,3]), 9)
  expect_equal(as.numeric(baf[6,9]), 2)
  expect_equal(sum(as.numeric(baf[6,-c(1,3,9)])), 0)
  # Gg
  expect_equal(as.numeric(baf[7,4]), 9)
  expect_equal(as.numeric(baf[7,10]), 2)
  expect_equal(sum(as.numeric(baf[7,-c(1,4,10)])), 0)
  # Tt
  expect_equal(as.numeric(baf[8,5]), 9)
  expect_equal(as.numeric(baf[8,11]), 2)
  expect_equal(sum(as.numeric(baf[8,-c(1,5,11)])), 0)
  
  # C^!
  expect_equal(as.numeric(baf[9,3]), 3)
  expect_equal(sum(as.numeric(baf[9,-c(1,3)])), 0)
})


##### ##### ##### ##### #####

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

