# Description: scMitoMut is a tool for detecting mitochondrial mutations in single-cell sequencing data

# Load libraries


#' scMitoMut
#'
#' @name scMitoMut
#' @docType package
#' @importFrom magrittr %>% %<>% extract extract2
#' @importFrom data.table data.table rbindlist .N :=
#' @importFrom plyr . count
#' @importFrom stringr str_glue str_match str_split_fixed str_c
#' @importFrom utils combn aregexec
#' @import stats
#' @importFrom Seurat DimPlot
#' @import rhdf5
#' @import readr
#' @import methods
#' @import Rcpp
#' @import ggplot2
#' @import zlibbioc
#' @useDynLib scMitoMut
NULL


