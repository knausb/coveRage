
library(covR)
context("Bedify")

ex_file <- system.file("extdata", "sc10_4k.mpileup.gz", package = "covR")
stats <- file_stats(ex_file, verbose=0)
x1 <- read_matrix(ex_file, nrows=stats[2], cols=1:stats[3], verbose=0)

# Create bed data
myBed <- matrix(ncol=4, nrow=5)
#myBed[,1] <- "sc12"
myBed[,1] <- "Supercontig_1.10"
myBed[,2] <- c(1, 1001, 1011, 9100, 10600)
myBed[,3] <- c(10, 1020, 1030, 9110, 10610)
myBed[,4] <- paste("gene", 1:5, sep="_")


# Without fill_missing.
myGenes <- bedify(myBed, x1, verbose=0)


test_that("myGenes contains the correct number of genes", {
  expect_equal(length(myGenes), 5)
})

test_that("myGenes contains the correct number of rows", {
  expect_equal(nrow(myGenes[[1]]), 10)
  expect_equal(nrow(myGenes[[2]]), 20)
  expect_equal(nrow(myGenes[[3]]), 20)
  expect_equal(nrow(myGenes[[4]]), 5)
  expect_equal(nrow(myGenes[[5]]), 8)
#  lapply(myGenes, nrow)
})



# With fill_missing.
myGenes <- bedify(myBed, x1[,c(1:3, 5, 6)], fill_missing=1)

#names(myGenes)


test_that("bedify return list names are correct", {
  expect_equal(sum(names(myGenes) == myBed[,4]), nrow(myBed))
})


test_that("myGenes contains the correct number of genes", {
  expect_equal(length(myGenes), 5)
})

test_that("myGenes contains the correct number of rows", {
  expect_equal(nrow(myGenes[[1]]), 10)
  expect_equal(nrow(myGenes[[2]]), 20)
  expect_equal(nrow(myGenes[[3]]), 20)
  expect_equal(nrow(myGenes[[4]]), 11)
  expect_equal(nrow(myGenes[[5]]), 11)
#  lapply(myGenes, nrow)
})

#test_that("myGenes contains NAs", {
#  expect_equal(sum(is.na(myGenes[[4]][1:6,3:54])), 312)
#})

#myGenes[[5]][9:11,]

#myGenes_bf <- lapply(myGenes, baf_stats)

#lapply(myGenes_bf, baf_summary)




