
library(coveRage)
context("baf_stats functions")

ex_file <- system.file("extdata", "sc10_4k.mpileup.gz", package = "covR")

#detach(package:covR, unload=T)

stats <- file_stats(ex_file, verbose=0)
x1 <- read_matrix(ex_file, nrows=stats[2], cols=1:stats[3], verbose=0)

#x2 <- x1[1:10,]
x2 <- x1

#RcppParallel::setThreadOptions(numThreads = 1)
#x2 <- baf_stats(calls=x1[,5], quals=x1[,6], ref=x1[,3], minq=0)


#test_that("baf_stats works",{
#  x2 <- baf_stats(calls=x1[,5], quals=x1[,6], ref=x1[,3], minq=0)
#  head(x2)
#})


#RcppParallel::setThreadOptions(numThreads = 2)

#head(x2)
#head(x2[,c(2, 3, 5, 6)])



#x3 <- baf_stats_st(x2[,c(1:3, 5, 6)], minq=0)
x3 <- baf_stats(x2[,c(1:3, 5, 6)], minq=0)

#head(x3)


#baf_plot(x3)
#baf_plot(x3, xlim=c(0,2000))
#baf_plot(x3, xlim=c(9001,9106))
#baf_plot(x3, xlim=c(9100,9110))

#count_df <- x3


myBed <- matrix(ncol=4, nrow=5)

myBed[,1] <- "sc12"

myBed[,2] <- c(1, 1001, 1011, 9100, 10600)
myBed[,3] <- c(10, 1020, 1030, 9110, 10607)
myBed[,4] <- paste("gene", 1:5, sep="_")
#myBed



myGenes <- bedify(myBed[1:5,], x3)
#myGenes

#baf_plot(myGenes[[1]])

#baf_plot(myGenes[[4]], na.rm=T)

#colnames(x2) <- paste("col", 1:ncol(x2), sep="_")

#myGenes2 <- bedify_sm(myBed, x2)
#myGenes2 <- bedify_sm(myBed, x2)



#lapply(myGenes2, function(x){x[,1:4]})

#x2[,2]

x4 <- cbind(10, x3)
colnames(x4)[1] <- "CHROM"

myGenes3 <- bedify_nm(myBed, x4)

#lapply(myGenes3, function(x){x[,1:13]})


