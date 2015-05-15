

#' @title baf_plot
#' @rdname baf_plot
#' @aliases baf_plot
#' 
#' @description BAlleleFreq plot.
#' 
#' @param count_df a data.frame containing info
#' @param ... arguments to be passed to methods
#' 
#' 
#' @references Laurie, Cathy C, Kimberly F Doheny, Daniel B Mirel, Elizabeth W Pugh, Laura J Bierut, Tushar Bhangale, Frederick Boehm, Neil E Caporaso, Marilyn C Cornelis, Howard J Edenberg and others.
#'   2010. Quality control and quality assurance in genotypic data for genome-wide association studies.
#'   Genetic Epidemiology 34(6): 591--602.
#' 
#' @export
baf_plot <- function(count_df, ...){
  
  count_df$POS <- as.numeric(count_df$POS)
  
  tot_count <- rowSums(count_df[,c('A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n')])
  
  fr_reads <- rowSums(count_df[,c('A', 'C', 'G', 'T', 'N')])
  fr_reads <- rbind(fr_reads, rowSums(count_df[,c('a', 'c', 'g', 't', 'n')]))
  fr_reads[2,] <- fr_reads[2,] * -1

  bxp_data <- c(fr_reads[1,], fr_reads[2,])
  bxp_cats <- rep(1:2, each=ncol(fr_reads))

  par(mar=c(0,0,0,0))
  par(oma=c(2,3,0.5,3))
  layout( mat=matrix(1:4, ncol=2, nrow=2, byrow=TRUE), widths= c(5, 1), heights=c(1, 3))

  #barplot(fr_reads, space=0, col="#1E90FF", border=NA, axes=FALSE)
  plot(1, type='n', xlim=range(count_df$POS), ylim=c(min(fr_reads[2,]), max(fr_reads[1,])), frame.plot=FALSE, axes=FALSE)
  points(count_df$POS, fr_reads[1,], type='h', col="#1E90FF")
  points(count_df$POS, fr_reads[2,], type='h', col="#00008B")
  axis(side=2, las=1)
  
  boxplot(bxp_data ~ bxp_cats, axes=FALSE, col=c("#1E90FF", "#00008B"))
  axis(side=4, las=1)

  plot(1, type="n", xlim=c(min(count_df$POS), max(count_df$POS)), ylim=c(0,1), frame.plot=FALSE, axes=FALSE)
  
  points(count_df$POS, unlist(count_df['A'])/tot_count, pch=20, col=3)
  points(count_df$POS, unlist(count_df['a'])/tot_count, pch=20, col=3)
  
  points(count_df$POS, unlist(count_df['C'])/tot_count, pch=20, col=5)
  points(count_df$POS, unlist(count_df['c'])/tot_count, pch=20, col=5)
  
  points(count_df$POS, unlist(count_df['G'])/tot_count, pch=20, col=1)
  points(count_df$POS, unlist(count_df['g'])/tot_count, pch=20, col=1)
  
  points(count_df$POS, unlist(count_df['T'])/tot_count, pch=20, col=2)
  points(count_df$POS, unlist(count_df['t'])/tot_count, pch=20, col=2)
  
  axis(side=2, at=c(0, 0.25, 0.33, 0.5, 0.66, 0.75, 1), las=1)
  axis(side=1)
  
  abline(h=c(0.25, 0.33, 0.5, 0.66, 0.75), lty=2)
  
  legend('topright', legend=c('A', 'C', 'G', 'T'), text.col=c(3, 5, 1, 2), bty="n", text.font=2)
  
  
  par(mar=c(5,4,4,2))
  par(mfrow=c(1,1))
  
}


