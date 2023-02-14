# Description: scMitoMut is a tool for detecting mitochondrial mutations in single-cell sequencing data

# Load libraries


#' scMitoMut
#'
#' @name scMitoMut
#' @docType package
#' @importFrom magrittr %>% %<>% extract extract2
#' @importFrom data.table data.table rbindlist .N :=
#' @importFrom plyr . count
#' @importFrom Ckmeans.1d.dp Ckmeans.1d.dp
#' @importFrom egg ggarrange
#' @importFrom stringr str_glue str_match str_split_fixed str_c
#' @importFrom utils combn aregexec
#' @import stats
#' @import methods
#' @import Rcpp
#' @import ggplot2
#' @import zlibbioc
#' @useDynLib CellBarcode
NULL


