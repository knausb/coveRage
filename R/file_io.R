


#'
#' @rdname file_io
#' @aliases file_stats
#' @title file_io
#' @name file_io
#' @import Rcpp
#'
#' @description Tools for rapid input of tabular data.
#'
#'
#' @details The function file_stats scans a file and returns the number of columns in the first row and the total number of rows.
#'
#'
#' @param filename name os file containing tabular data
#' @param sep character which delimits columns
#' @param skip number of rows to skip before reading data
#' 
#' @export
file_stats <- function(filename, sep = '\t', skip = 0L){
    .Call('covR_file_stats', PACKAGE = 'covR', filename, sep, skip)
}


