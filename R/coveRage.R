#' coveRage: A package for exploring coverage in high throughput sequencing projects.
#'
#' This package provides functions to process mpileup files
#' created by samtools mpileup.  This information can be
#' used to explore copy number variants to help determine
#' whether they are due to biological or technical issues.
#' These data can then be used to censor genome analyses
#' to focus on regions of base ploid. 
#' 
#' 
#' @section File input and output functions:
#' file_stats
#' 
#' read_matrix
#' 
#' write_matrix
#'
#' @section Mpileup processing functions:
#' bedify
#' 
#' baf_stats
#' 
#' baf_plot
#' 
#' 
#' @section Vignettes:
#' Package vignettes can be listed with:
#' 
#' vignette(package='coveRage')
#' 
#' @seealso
#' \href{http://www.htslib.org/doc/#manuals}{samtools documentation} (and others).
#'
#'
#' @docType package
#' @name coveRage
#' @useDynLib coveRage
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' 
#' 
#' 
NULL
#> NULL