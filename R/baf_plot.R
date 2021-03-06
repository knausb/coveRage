

#' @title baf_plot
#' @rdname baf_plot
#' @aliases baf_plot
#' 
#' @description B allele frequency plot.
#' 
#' @param counts a numeric matrix containing count data
#' @param alpha opacity (0-255)
#' @param title for the plot.
#' @param na.rm Logical, should sites which contain NAs (no reported coverage) be removed
#' @param aline coordinates at which to draw tick marks on the y axis
#' @param hline coordinates at which to draw horizontal lines where NULL omits lines
#' @param halpha opacity (0-255) for horizontal lines
#' @param vline coordinates at which to draw vertical lines where NULL omits lines
#' @param valpha opacity (0-255) for vertical lines
#' @param ... arguments to be passed to methods
#' 
#' @details
#' Plots coverage and allele frequencies for a chromosomal region.
#' The upper pane presents per site sequencing depth while the lower pane presents the frequency of all four alleles.
#' These plots can be used to identify technical error (i.e., read alignment issues) and/or biological phenomena such as copy number variation.
#' 
#' The matrix of count data can be generated from mpileup data with the function \code{\link{baf_stats}}
#' 
#' @seealso \code{\link{baf_stats}}
#' 
#' @references Laurie, Cathy C, Kimberly F Doheny, Daniel B Mirel, Elizabeth W Pugh, Laura J Bierut, Tushar Bhangale, Frederick Boehm, Neil E Caporaso, Marilyn C Cornelis, Howard J Edenberg and others.
#'   2010. Quality control and quality assurance in genotypic data for genome-wide association studies.
#'   Genetic Epidemiology 34(6): 591--602.
#' 
#' @export
baf_plot <- function(counts, 
                     alpha=255, 
                     title="Locus", 
                     na.rm=FALSE,
                     aline=c(0, 0.25, 0.33, 0.5, 0.66, 0.75, 1),
                     hline=c(0.25, 0.33, 0.5, 0.66, 0.75),
                     halpha=255,
                     vline=NULL,
                     valpha=255,
                     ...){
  
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

  graphics::par(mar=c(0,0,0,0))
#  graphics::par(oma=c(3,4,0.5,3))
  graphics::par(oma=c(3,4,1.5,3))
#  graphics::layout( mat=matrix(1:4, ncol=2, nrow=2, byrow=TRUE), widths= c(8, 1), heights=c(1, 3))
  graphics::layout( mat=matrix(c(1,2,4,3), ncol=2, nrow=2, byrow=TRUE), widths= c(8, 1), heights=c(1, 3))
  
  # Depth plot.
  #barplot(fr_reads, space=0, col="#1E90FF", border=NA, axes=FALSE)
  graphics::plot(range(counts[,'POS']), c( max(fr_reads[1,]), min(fr_reads[2,]) ), type='n', frame.plot=FALSE, axes=FALSE, ylab="", xlab="", ...)
#  plot(1, type='n', xlim=range(counts$POS), ylim=c(min(fr_reads[2,]), max(fr_reads[1,])), frame.plot=FALSE, axes=FALSE)
  graphics::abline(h=0)
  graphics::points(counts[,'POS'], fr_reads[1,], type='h', col="#1E90FF")
  graphics::points(counts[,'POS'], fr_reads[2,], type='h', col="#00008B")
  graphics::axis(side=2, las=1)
  graphics::mtext(text="F/R coverage", side=2, line=3)

  # Marginal boxplots.
  graphics::boxplot(bxp_data ~ bxp_cats, axes=FALSE, col=c("#1E90FF", "#00008B"))
  graphics::axis(side=4, las=1)

  # Legend.
  graphics::plot(1, axes=FALSE, type='n', xlab="", ylab="")
  graphics:: legend('center', legend=c('A', 'C', 'G', 'T'), text.col=c(3, 5, 1, 2), bty="n", text.font=2, cex=2, xjust=0.5)
  title(main=title, outer=T)
  
  # BAF plot.
  graphics::plot(range(counts[,'POS']), c(0,1), type="n", frame.plot=FALSE, axes=FALSE, ylab="", xlab="", ...)
#  plot(1, type="n", xlim=c(min(counts$POS), max(counts$POS)), ylim=c(0,1), frame.plot=FALSE, axes=FALSE)

  # Sum over forward and reverse reads.
  graphics::points(counts[,'POS'], rowSums(counts[,c('A','a')])/tot_count, pch=20, col=grDevices::rgb(0, 205, 0, alpha, maxColorValue=255))
  graphics::points(counts[,'POS'], rowSums(counts[,c('C','c')])/tot_count, pch=20, col=grDevices::rgb(0, 255, 255, alpha, maxColorValue=255))
  graphics::points(counts[,'POS'], rowSums(counts[,c('G','g')])/tot_count, pch=20, col=grDevices::rgb(0, 0, 0, alpha, maxColorValue=255))
  graphics::points(counts[,'POS'], rowSums(counts[,c('T','t')])/tot_count, pch=20, col=grDevices::rgb(255, 0, 0, alpha, maxColorValue=255))

  # Use only forward reads.
#  graphics::points(counts[,'POS'], counts[ ,c('A') ]/tot_count, pch=20, col=grDevices::rgb(0, 205, 0, alpha, maxColorValue=255))
#  graphics::points(counts[,'POS'], counts[ ,c('C') ]/tot_count, pch=20, col=grDevices::rgb(0, 255, 255, alpha, maxColorValue=255))
#  graphics::points(counts[,'POS'], counts[ ,c('G') ]/tot_count, pch=20, col=grDevices::rgb(0, 0, 0, alpha, maxColorValue=255))
#  graphics::points(counts[,'POS'], counts[ ,c('T') ]/tot_count, pch=20, col=grDevices::rgb(255, 0, 0, alpha, maxColorValue=255))


  graphics::axis(side=2, at=aline, las=1)
  graphics::axis(side=1)
  graphics::abline(h=hline, lty=2, col=grDevices::rgb(0, 0, 0, halpha, maxColorValue = 255))
  graphics::abline(v=vline, lty=2, col=grDevices::rgb(0, 0, 0, valpha, maxColorValue = 255))
  
  graphics::mtext(text="Allele frequency", side=2, line=3)
  graphics::mtext(text="Position", side=1, line=2)


  graphics::par(mar=c(5,4,4,2))
  graphics::par(mfrow=c(1,1))
  
}


