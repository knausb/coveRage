

#' @title baf_plot
#' @rdname baf_plot
#' @aliases baf_plot
#' 
#' @description BAlleleFreq plot.
#' 
#' @param counts a numeric matrix containing count data
#' @param alpha opacity (0-255)
#' @param na.rm Logical, should sites which contain NAs be removed
#' @param ... arguments to be passed to methods
#' 
#' 
#' @references Laurie, Cathy C, Kimberly F Doheny, Daniel B Mirel, Elizabeth W Pugh, Laura J Bierut, Tushar Bhangale, Frederick Boehm, Neil E Caporaso, Marilyn C Cornelis, Howard J Edenberg and others.
#'   2010. Quality control and quality assurance in genotypic data for genome-wide association studies.
#'   Genetic Epidemiology 34(6): 591--602.
#' 
#' @export
baf_plot <- function(counts, alpha=255, na.rm=FALSE, ...){
  
#  if(class(counts$POS) != "numeric"){
#    counts$POS <- as.numeric(counts$POS)
#  }
  
#  counts$POS <- as.numeric(counts$POS)
  
  if(na.rm==TRUE){
    counts <- counts[!colSums(apply(counts, MARGIN=1, is.na))>0,]
  }
  
  tot_count <- rowSums(counts[,c('A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n')])
  
  fr_reads <- rowSums(counts[,c('A', 'C', 'G', 'T', 'N')])
  fr_reads <- rbind(fr_reads, rowSums(counts[,c('a', 'c', 'g', 't', 'n')]))
  fr_reads[2,] <- fr_reads[2,] * -1

  bxp_data <- c(fr_reads[1,], fr_reads[2,])
  bxp_cats <- rep(1:2, each=ncol(fr_reads))

  par(mar=c(0,0,0,0))
  par(oma=c(3,4,0.5,3))
  layout( mat=matrix(1:4, ncol=2, nrow=2, byrow=TRUE), widths= c(8, 1), heights=c(1, 3))

  #barplot(fr_reads, space=0, col="#1E90FF", border=NA, axes=FALSE)
  plot(range(counts[,'POS']), c( max(fr_reads[1,]), min(fr_reads[2,]) ), type='n', frame.plot=FALSE, axes=FALSE, ylab="", xlab="", ...)
#  plot(1, type='n', xlim=range(counts$POS), ylim=c(min(fr_reads[2,]), max(fr_reads[1,])), frame.plot=FALSE, axes=FALSE)
  abline(h=0)
  points(counts[,'POS'], fr_reads[1,], type='h', col="#1E90FF")
  points(counts[,'POS'], fr_reads[2,], type='h', col="#00008B")
  axis(side=2, las=1)
  mtext(text="F/R coverage", side=2, line=3)

  boxplot(bxp_data ~ bxp_cats, axes=FALSE, col=c("#1E90FF", "#00008B"))
  axis(side=4, las=1)


  plot(range(counts[,'POS']), c(0,1), type="n", frame.plot=FALSE, axes=FALSE, ylab="", xlab="", ...)
#  plot(1, type="n", xlim=c(min(counts$POS), max(counts$POS)), ylim=c(0,1), frame.plot=FALSE, axes=FALSE)

  points(counts[,'POS'], rowSums(counts[,c('A','a')])/tot_count, pch=20, col=rgb(0, 205, 0, alpha, maxColorValue=255))
  points(counts[,'POS'], rowSums(counts[,c('C','c')])/tot_count, pch=20, col=rgb(0, 255, 255, alpha, maxColorValue=255))
  points(counts[,'POS'], rowSums(counts[,c('G','g')])/tot_count, pch=20, col=rgb(0, 0, 0, alpha, maxColorValue=255))
  points(counts[,'POS'], rowSums(counts[,c('T','t')])/tot_count, pch=20, col=rgb(255, 0, 0, alpha, maxColorValue=255))

  axis(side=2, at=c(0, 0.25, 0.33, 0.5, 0.66, 0.75, 1), las=1)
  axis(side=1)
  
  abline(h=c(0.25, 0.33, 0.5, 0.66, 0.75), lty=2)
  mtext(text="Allele frequency", side=2, line=3)
  mtext(text="Position", side=1, line=2)

  plot(1, axes=FALSE, type='n', xlab="", ylab="")

  legend('center', legend=c('A', 'C', 'G', 'T'), text.col=c(3, 5, 1, 2), bty="n", text.font=2, cex=2, xjust=0.5)

  par(mar=c(5,4,4,2))
  par(mfrow=c(1,1))
  
}


