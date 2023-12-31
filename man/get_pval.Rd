% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{get_pval}
\alias{get_pval}
\title{Get p-value list for single locus}
\usage{
get_pval(mtmutObj, loc, model = "bb", method = "fdr")
}
\arguments{
\item{mtmutObj}{a mtmutObj object.}

\item{loc}{a string of the locus.}

\item{model}{a string of the model for mutation calling, it can be "bb", "bm" or "bi" which stands for beta binomial, binomial mixture and binomial model respectively.}

\item{method}{a string of the method for p-value adjustment, refer to \code{\link[stats]{p.adjust}}.}
}
\value{
a vector of p-value for each cell.
}
\description{
This function returns the p-value list for a single locus.
}
\examples{
## Use the example data
f <- system.file("extdata", "mini_dataset.tsv.gz", package = "scMitoMut")
## Create a temporary h5 file
## In real case, we keep the h5 in project folder for future use
f_h5_tmp <- tempfile(fileext = ".h5")

## Load the data with parse_table function
f_h5 <- parse_table(f, sep = "\t", h5_file = f_h5_tmp)
## open the h5 file and create a mtmutObj object
x <- open_h5_file(f_h5)
run_model_fit(x)
get_pval(x, "chrM.1000", "bb", "fdr")
get_pval(x, "chrM.1000", "bm", "fdr")
get_pval(x, "chrM.1000", "bi", "fdr")
}
