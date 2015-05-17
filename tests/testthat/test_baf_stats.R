

context("baf_stats functions")

ex_file <- system.file("extdata", "sc10_4k.mpileup.gz", package = "covR")

#detach(package:covR, unload=T)
library(covR)
stats <- file_stats(ex_file, verbose=0)
x1 <- read_matrix(ex_file, nrows=stats[2], ncols=stats[3], verbose=0)

#x2 <- x1[1:10,]
x2 <- x1

#RcppParallel::setThreadOptions(numThreads = 1)
#x2 <- baf_stats(calls=x1[,5], quals=x1[,6], ref=x1[,3], minq=0)


#test_that("baf_stats works",{
#  x2 <- baf_stats(calls=x1[,5], quals=x1[,6], ref=x1[,3], minq=0)
#  head(x2)
#})


#RcppParallel::setThreadOptions(numThreads = 2)



x3 <- baf_stats_st(calls=x2[,5], quals=x2[,6], ref=x2[,3], minq=0)
#x3 <- baf_stats_st(calls=x2[,5], quals=x2[,6], ref=x2[,3], minq=30)

#x3
#x2[,c(3,5)]

#x3 <- cbind(as.character(x2[,1]), as.numeric(x2[,2]), as.data.frame(x3))


#names(x3)[1] <- "CHROM"
#names(x3)[2] <- "POS"
#x3$POS <- as.numeric(x3$POS)

x3 <- cbind(as.numeric(x2[,2]), x3)
colnames(x3)[1] <- "POS"

head(x3)


baf_plot(x3)
baf_plot(x3, xlim=c(0,2000))

#count_df <- x3


myBed <- matrix(ncol=4, nrow=3)

myBed[,1] <- "sc12"
myBed[,2] <- c(1001, 1011, 9001)
myBed[,3] <- c(1020, 1030, 9010)
myBed[,4] <- paste("gene", 1:3, sep="_")
myBed


myGenes <- bedify(myBed, x3)
myGenes
