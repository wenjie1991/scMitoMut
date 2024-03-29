######################################################################
#                          Internal function                          #
#######################################################################

## Parse the output of mtGATK into a data.table.
read_mgatk <- function(mgatk_output_dir, prefix) {
    mgatk_output_dir <- paste0(mgatk_output_dir)
    base_reads_files <- dir(mgatk_output_dir, str_glue("{prefix}.[ACTG].txt.gz"), full.names = TRUE)

    message("Allele count files are:")
    message(base_reads_files)

    names(base_reads_files) <- str_split_fixed(base_reads_files %>% basename(), "\\.", 4)[, 2]
    ref_file <- dir(mgatk_output_dir, "ref", full.names = TRUE)
    coverage_file <- dir(mgatk_output_dir, "coverage", full.names = TRUE)

    ## base-wise sequence depth
    base_depth_d <- lapply(seq_along(base_reads_files), function(i) {
        f <- base_reads_files[i]
        d <- data.table::fread(f)
        names(d) <- c("loc", "cell_barcode", "fwd_depth", "rev_depth")
        d$alt <- names(base_reads_files)[i]
        d
}) %>% rbindlist()

    ## reference sequence
    ref_d <- data.table::fread(ref_file)
    names(ref_d) <- c("loc", "ref")

    ## base-wise coverage
    if (length(coverage_file) == 0) {
        coverage_d <- base_depth_d[, .(coverage = sum(fwd_depth + rev_depth, na.rm = TRUE)), by = c("loc", "cell_barcode")]
    } else {
        coverage_d <- data.table::fread(coverage_file)
        names(coverage_d) <- c("loc", "cell_barcode", "coverage")
    }

    ## merge base-wise depth, ref seq, coverage
    merge_d <- merge(base_depth_d, coverage_d, by = c("loc", "cell_barcode"))
    merge_d <- merge(merge_d, ref_d, by = "loc")
    # rm(base_depth_d)
    # rm(coverage_d)
    # rm(ref_d)
    # gc()
    merge_d
}

## This function sorts the matrix for better visualization of mutual exclusivity across genes
## Adaptoed from: https://gist.github.com/armish/564a65ab874a770e2c26
memoSort <- function(M, m = NULL) {
    geneOrder <- sort(rowSums(M, na.rm = TRUE), decreasing = TRUE, index.return = TRUE)$ix
    scoreCol <- function(x) {
        score <- 0
        for (i in seq_along(x)) {
            if (x[i] & !is.na(x[i])) {
                score <- score + 2^(length(x) - i)
            }
        }
        return(score)
    }
    scores <- apply(M[geneOrder, ], 2, scoreCol)
    sampleOrder <- sort(scores, decreasing = TRUE, index.return = TRUE)$ix
    if (is.null(m)) {
        return(M[geneOrder, sampleOrder])
    } else {
        return(m[geneOrder, sampleOrder])
    }
}

## Read one locus
read_locus <- function(mtmutObj, loc, maj_base = NULL) {
    d_sub <- data.table((mtmutObj$mut_table & loc)[])[cell_barcode %in% mtmutObj$cell_selected]

    ## maj_base is the base with max frequency
    if (is.null(maj_base)) {
        # base_f <- table(d_sub$alt)
        # maj_base <- names(base_f)[which.max(base_f)]
        maj_base <- d_sub[,
            .(af = median((fwd_depth + rev_depth) / coverage)),
            by = alt
            ][which.max(af), alt]
    }

    d_coverage <- unique(d_sub[, .(cell_barcode, coverage)])

    d_select_maj_base <- d_sub[alt == maj_base, .(cell_barcode, alt_depth = fwd_depth + rev_depth, fwd_depth, rev_depth)]
    d_select_maj_base <- merge(d_select_maj_base, d_coverage, by = "cell_barcode", all = TRUE)
    d_select_maj_base[is.na(d_select_maj_base)] <- 0

    ## TODO: do we need this?
    # d_select_maj_base = d_select_maj_base[alt_depth / coverage > 0.8]

    d_select_maj_base
}


#######################################################################
#                    End internal function region                     #
#######################################################################


#' Load allele count table
#'
#' This function loads the allele count table and save it to a H5 file.
#'
#' @param file a string of the allele count table file directory.
#' @param h5_file a string of the output h5 file directory.
#' @param ... other parameters passed to \code{\link[data.table]{fread}}.
#' @return a string of the output h5 file directory.
#' @details The allele count table should be a data.table with the following columns:
#' \describe{
#'  \item{loc}{a string of the locus}
#'  \item{cell_barcode}{a string of the cell barcode.}
#'  \item{fwd_depth}{a integer of the forward depth.}
#'  \item{rev_depth}{a integer of the reverse depth.}
#'  \item{alt}{a string of the alternative base.}
#'  \item{ref}{a string of the reference base.}
#'  \item{coverage}{a integer of the coverage.}
#' }
#'
#' @examples
#' ## Use the example data
#' f <- system.file("extdata", "mini_dataset.tsv.gz", package = "scMitoMut")
#'
#' ## Create a temporary h5 file
#' ## In real case, we keep the h5 in project folder for future use
#' f_h5_tmp <- tempfile(fileext = ".h5")
#'
#' ## Load the data with parse_table function
#' f_h5 <- parse_table(f, sep = "\t", h5_file = f_h5_tmp)
#' f_h5
#'
#' ## open the h5 file and create a mtmutObj object
#' x <- open_h5_file(f_h5)
#' x
#' @export
parse_table <- function(file, h5_file = "mut.h5", ...) {

    ## check file
    if (!file.exists(file)) {
        stop("File not exists!")
    }

    merge_d <- data.table::fread(file, ...)

    if (!all(c("loc", "cell_barcode", "fwd_depth", "rev_depth", "alt", "ref", "coverage") %in% names(merge_d))) {
        stop("The allele count table should have the following columns: loc, cell_barcode, fwd_depth, rev_depth, alt, ref, coverage")
    }



    ##############################
    ## save to h5 file
    if (file.exists(h5_file)) {
        warning("\nH5 file exists, remove it.\n")
        file.remove(h5_file)
    }

    H5Fcreate(h5_file)
    h5f <- H5Fopen(h5_file)

    ## mut_table
    h5g <- H5Gcreate(h5loc = h5f, name = "mut_table")
    plyr::d_ply(merge_d, "loc", function(x) {
        d_name <- paste0("chrM.", x$loc[1])
        x$loc <- NULL
        h5write(x, h5g, d_name)
})

    ## cell list
    cell_list <- merge_d$cell_barcode %>% unique()
    h5write(cell_list, h5f, "cell_list")
    h5write(cell_list, h5f, "cell_selected")

    ## loc list
    loc_list <- merge_d$loc %>%
        unique() %>%
        paste0("chrM.", .)
    h5write(loc_list, h5f, "loc_list")
    h5write(loc_list, h5f, "loc_selected")

    ## TODO: Do we need to close the H5 file?
    #     H5Fclose(h5f)
    #     H5Gclose(h5g)
    #     h5closeAll()

    ## Remove the big data
    rm("merge_d")
    gc()

    h5_file
}


#' Load mtGATK output
#'
#' This function loads the mtGATK output and save it to a H5 file.
#'
#' @param dir a string of the mtGATK output \code{final} fold.
#' @param prefix a string of the prefix of the mtGATK output directory.
#' @param h5_file a string of the output h5 file directory.
#' @return a string of the output h5 file directory.
#'
#' @examples
#' ## Use the allele count table data
#' f <- system.file("extdata", "mini_dataset.tsv.gz", package = "scMitoMut")
#' ## Create a temporary h5 file
#' ## In real case, we keep the h5 in project folder for future use
#' f_h5_tmp <- tempfile(fileext = ".h5")
#' ## Load the data with parse_table function
#' f_h5 <- parse_table(f, sep = "\t", h5_file = f_h5_tmp)
#' f_h5
#' ## Use the mgatk output
#' f <- system.file("extdata", "mini_mgatk_out", package = "scMitoMut")
#' f_h5_tmp <- tempfile(fileext = ".h5")
#' f_h5 <- parse_mgatk(paste0(f, "/final/"), prefix = "sample", h5_file = f_h5_tmp)
#' f_h5
#' x <- open_h5_file(f_h5)
#' x
#' ##
#' @export
parse_mgatk <- function(dir, prefix, h5_file = "mut.h5") {

    ## check dir 
    if (!dir.exists(dir)) {
        stop("Directory ", dir, " not exists!")
    }

    ##############################
    ## Read in data
    merge_d <- read_mgatk(mgatk_output_dir = dir, prefix = prefix)

    ##############################
    ## save to h5 file
    if (file.exists(h5_file)) {
        warning("\nH5 file exists, remove it.\n")
        file.remove(h5_file)
    }

    H5Fcreate(h5_file)
    h5f <- H5Fopen(h5_file)

    ## mut_table
    h5g <- H5Gcreate(h5loc = h5f, name = "mut_table")
    plyr::d_ply(merge_d, "loc", function(x) {
        d_name <- paste0("chrM.", x$loc[1])
        x$loc <- NULL
        h5write(x, h5g, d_name)
})

    ## cell list
    cell_list <- merge_d$cell_barcode %>% unique()
    h5write(cell_list, h5f, "cell_list")
    h5write(cell_list, h5f, "cell_selected")

    ## loc list
    loc_list <- merge_d$loc %>%
        unique() %>%
        paste0("chrM.", .)
    h5write(loc_list, h5f, "loc_list")
    h5write(loc_list, h5f, "loc_selected")

    ## TODO: Do we need to close the H5 file?
    #     H5Fclose(h5f)
    #     H5Gclose(h5g)
    #     h5closeAll()

    ## Remove the big data
    rm("merge_d")
    gc()

    h5_file
}

#' Open H5 file
#'
#' This function opens the H5 file and create a mtmutObj object.
#'
#' @param h5_file a string of the h5 file directory
#' @return a mtmutObj object
#' @details The mtmutObj object is a S3 class for handling mitochondrial mutation data.
#' It contains the following elements:
#' \describe{
#'  \item{file}{a string of the h5 file directory.}
#'  \item{h5f}{H5 file handle.}
#'  \item{mut_table}{allele count table H5 group handle.}
#'  \item{loc_list}{list of available loci.}
#'  \item{loc_selected}{selected loci, the default is all loci.}
#'  \item{cell_list}{list of available cell ids.}
#'  \item{cell_selected}{selected cell ids, the default is all cells.}
#'  \item{loc_pass}{loci passed the filter, the default is NULL}
#'  \item{loc_filter}{filter parameters.}
#'  \item{loc_filter$min_cell}{a integer of the minimum number of cells with mutation, the default is 1.}
#'  \item{loc_filter$model}{a string of the model for mutation calling, it can be "bb", "bm" or "bi", the default is "bb".}
#'  \item{loc_filter$p_threshold}{a numeric of the p-value threshold, the default is 0.05.}
#'  \item{loc_filter$p_adj_method}{a string of the method for p-value adjustment, refer to \code{\link[stats]{p.adjust}}, the default is "fdr".}
#'  }
#'
#' @examples
#' ## Use the example data
#' f <- system.file("extdata", "mini_dataset.tsv.gz", package = "scMitoMut")
#' ## Create a temporary h5 file
#' ## In real case, we keep the h5 in project folder for future use
#' f_h5_tmp <- tempfile(fileext = ".h5")
#' ## Load the data with parse_table function
#' f_h5 <- parse_table(f, sep = "\t", h5_file = f_h5_tmp)
#' f_h5
#' ## open the h5 file and create a mtmutObj object
#' x <- open_h5_file(f_h5)
#' x
#' @export
open_h5_file <- function(h5_file) {

    ## check file
    if (!file.exists(h5_file)) {
        stop("File ", h5_file, " not exists!")
    }
    
    h5f <- H5Fopen(h5_file)
    h5g <- h5f & "mut_table"
    loc_list <- h5f$loc_list
    loc_selected <- h5f$loc_selected
    cell_list <- h5f$cell_list
    cell_selected <- h5f$cell_selected
    ## we call the output as mtmutObj
    mtmutObj <- list(
        file = h5_file,
        h5f = h5f,
        mut_table = h5g,
        loc_list = loc_list,
        loc_selected = loc_selected,
        cell_list = cell_list,
        cell_selected = cell_selected,
        loc_pass = NULL,
        loc_filter = list(
            min_cell = 1,
            model = "bb",
            p_threshold = 0.05,
            alt_count_threshold = 0,
            p_adj_method = "fdr"
        )
    )
    class(mtmutObj) <- "mtmutObj"
    mtmutObj
}

#' @export
#' @rdname print.mtmutObj
format.mtmutObj <- function(x, ...) {
    cat("mtmutObj object\n")
    cat("-------------------------------------------------\n")
    cat("h5 file: ")
    cat(x$file)
    cat("\nAvailable loci: ")
    cat(length(x$loc_list))
    cat("\nSelected loci: ")
    cat(length(x$loc_selected))
    cat("\nAvailable cells: ")
    cat(length(x$cell_list))
    cat("\nSelected cells: ")
    cat(length(x$cell_selected))
    cat("\nLoci passed the filter: ")
    cat(length(x$loc_pass))
    cat("\nfilter parameters: ")
    cat("\n")
    cat("\t", "min_cell: ", x$loc_filter$min_cell, "\n", sep = "")
    cat("\t", "model: ", x$loc_filter$model, "\n", sep = "")
    cat("\t", "p_threshold: ", x$loc_filter$p_threshold, "\n", sep = "")
    cat("\t", "alt_count_threshold: ", x$loc_filter$alt_count_threshold, "\n", sep = "")
    cat("\t", "p_adj_method: ", x$loc_filter$p_adj_method, "\n", sep = "")
}

#' Print mtmutObj object
#'
#' The print method for mtmutObj object.
#'
#' @param x a mtmutObj object.
#' @param ... other parameters passed to \code{\link[base]{format}} or \code{\link[base]{print}}.
#' @return a string
#' @export
#' @examples
#' ## Use the example data
#' f <- system.file("extdata", "mini_dataset.tsv.gz", package = "scMitoMut")
#' ## Create a temporary h5 file
#' ## In real case, we keep the h5 in project folder for future use
#' f_h5_tmp <- tempfile(fileext = ".h5")
#' ## Load the data with parse_table function
#' f_h5 <- parse_table(f, sep = "\t", h5_file = f_h5_tmp)
#' f_h5
#' ## open the h5 file and create a mtmutObj object
#' x <- open_h5_file(f_h5)
#' x
#' print(x)
print.mtmutObj <- function(x, ...) {
    format(x, ...)
}

#' @export
#' @rdname print.mtmutObj
is.mtmutObj <- function(x) inherits(x, "mtmutObj")


#' Subset cell and loci
#'
#' Functions to subset cell and loci for fitting models and plotting.
#'
#' @param mtmutObj a mtmutObj object
#' @param cell_list a list of cell barcodes
#' @param loc_list a list of loci
#' @return a mtmutObj object with cell and loci selected
#' @examples
#' ## Use the example data
#' f <- system.file("extdata", "mini_dataset.tsv.gz", package = "scMitoMut")
#'
#' ## Create a temporary h5 file
#' ## In real case, we keep the h5 in project folder for future use
#' f_h5_tmp <- tempfile(fileext = ".h5")
#'
#' ## Load the data with parse_table function
#' f_h5 <- parse_table(f, sep = "\t", h5_file = f_h5_tmp)
#'
#' ## open the h5 file and create a mtmutObj object
#' x <- open_h5_file(f_h5)
#' x
#' ## subset cell and loci
#' x <- subset_cell(x, x$cell_list[1:10])
#' x <- subset_loc(x, x$loc_list[1:10])
#' x
#' @export
subset_cell <- function(mtmutObj, cell_list) {
    if (!is(mtmutObj, "mtmutObj")) {
        stop("mtmutObj should be a mtmutObj object")
    }

    if ("cell_selected" %in% h5ls(mtmutObj$h5f, recursive = FALSE)$name) {
        h5delete(mtmutObj$h5f, "cell_selected")
    }
    if (length(unique(cell_list)) != length(cell_list)) {
        stop("cell_list should be unique")
    }
    h5write(cell_list, mtmutObj$h5f, "cell_selected")
    mtmutObj$cell_selected <- cell_list
    mtmutObj
}

#' @rdname subset_cell
#' @export
subset_loc <- function(mtmutObj, loc_list) {
    if (!is(mtmutObj, "mtmutObj")) {
        stop("mtmutObj should be a mtmutObj object")
    }

    if ("loc_selected" %in% h5ls(mtmutObj$h5f, recursive = FALSE)$name) {
        h5delete(mtmutObj$h5f, "loc_selected")
    }

    if (length(unique(loc_list)) != length(loc_list)) {
        stop("loc_list should be unique")
    }

    h5write(loc_list, mtmutObj$h5f, "loc_selected")
    mtmutObj$loc_selected <- loc_list
    mtmutObj
}

#' Get p-value list for single locus
#'
#' This function returns the p-value list for a single locus.
#'
#' @param mtmutObj a mtmutObj object.
#' @param loc a string of the locus.
#' @param model a string of the model for mutation calling, it can be "bb", "bm" or "bi" which stands for beta binomial, binomial mixture and binomial model respectively.
#' @param method a string of the method for p-value adjustment, refer to \code{\link[stats]{p.adjust}}.
#' @return a vector of p-value for each cell.
#' @examples
#' ## Use the example data
#' f <- system.file("extdata", "mini_dataset.tsv.gz", package = "scMitoMut")
#' ## Create a temporary h5 file
#' ## In real case, we keep the h5 in project folder for future use
#' f_h5_tmp <- tempfile(fileext = ".h5")
#'
#' ## Load the data with parse_table function
#' f_h5 <- parse_table(f, sep = "\t", h5_file = f_h5_tmp)
#' ## open the h5 file and create a mtmutObj object
#' x <- open_h5_file(f_h5)
#' run_model_fit(x)
#' get_pval(x, "chrM.1000", "bb", "fdr")
#' get_pval(x, "chrM.1000", "bm", "fdr")
#' get_pval(x, "chrM.1000", "bi", "fdr")
#' @export
get_pval <- function(mtmutObj, loc, model = "bb", method = "fdr") {
    if (!is(mtmutObj, "mtmutObj")) {
        stop("mtmutObj should be a mtmutObj object")
    }

    if (model == "bb") {
        pval_item <- "bb_pval"
    } else if (model == "bm") {
        pval_item <- "bm_pval"
    } else if (model == "bi") {
        pval_item <- "bi_pval"
    } else {
        stop("model should be bb, bm or bi")
    }
    p.adjust((mtmutObj$h5f & "pval" & loc & pval_item)[], method)
}

#' Filter mutations
#'
#' This function filters the mutations based on the mutation calling model and parameters. The loci passed the filter will be saved in the h5 file, together with the filter parameters.
#'
#' @param mtmutObj a mtmutObj object.
#' @param min_cell a integer of the minimum number of cells with mutation, the default is 5.
#' @param model a string of the model for mutation calling, it can be "bb", "bm" or "bi" which stands for beta binomial, binomial mixture and binomial model respectively, the default is "bb".
#' @param p_threshold a numeric of the p-value threshold, the default is 0.05.
#' @param p_adj_method a string of the method for p-value adjustment, .
#'   refer to \code{\link[stats]{p.adjust}}. The default is "fdr".
#' @param alt_count_threshold a integer of the minimum number of alternative base count, the default is 0.
#' @return a mtmutObj object with loc_pass and loc_filter updated.
#' @examples
#' ## Use the example data
#' f <- system.file("extdata", "mini_dataset.tsv.gz", package = "scMitoMut")
#'
#' ## Create a temporary h5 file
#' ## In real case, we keep the h5 in project folder for future use
#' f_h5_tmp <- tempfile(fileext = ".h5")
#'
#' ## Load the data with parse_table function
#' f_h5 <- parse_table(f, sep = "\t", h5_file = f_h5_tmp)
#'
#' ## open the h5 file and create a mtmutObj object
#' x <- open_h5_file(f_h5)
#' run_model_fit(x)
#' x <- filter_loc(x, min_cell = 5, model = "bb", p_threshold = 0.05, p_adj_method = "fdr")
#' x
#' @export
filter_loc <- function(mtmutObj, min_cell = 5, model = "bb", p_threshold = 0.05, alt_count_threshold = 0, p_adj_method = "fdr") {

    if (!is(mtmutObj, "mtmutObj")) {
        stop("mtmutObj should be a mtmutObj object")
    }

    if (min_cell < 0) {
        stop("min_cell should be >= 1")
    }

    if (model != "bb" & model != "bm" & model != "bi") {
        stop("model should be bb, bm or bi")
    }

    if (p_threshold < 0 | p_threshold > 1) {
        stop("p_threshold should be in [0, 1]")
    }

    if (alt_count_threshold < 0) {
        stop("alt_count_threshold should be >= 0")
    }

    ## Get the parameters from mtmutObj if they are NULL
    min_cell <- ifelse(is.null(min_cell), mtmutObj$loc_filter$min_cell, min_cell)
    model <- ifelse(is.null(model), mtmutObj$loc_filter$model, model)
    p_threshold <- ifelse(is.null(p_threshold), mtmutObj$loc_filter$p_threshold, p_threshold)
    alt_count_threshold <- ifelse(is.null(alt_count_threshold), mtmutObj$loc_filter$alt_count_threshold, alt_count_threshold)
    p_adj_method <- ifelse(is.null(p_adj_method), mtmutObj$loc_filter$p_adj_method, p_adj_method)

    loc_list <- mtmutObj$loc_selected
    res <- parallel::mclapply(loc_list, function(xi) {
        pval <- get_pval(mtmutObj, xi, model = model, method = p_adj_method)
        data.frame(loc = xi, mut_cell_n = sum(pval <= p_threshold, na.rm = TRUE))
    }) %>% rbindlist()
    res <- res[mut_cell_n >= min_cell]

    if (alt_count_threshold > 0) {
        res_2 <- parallel::mclapply(res$loc, function(xi) {
            d <- read_locus(mtmutObj, xi)
            pval <- get_pval(mtmutObj, xi, model = model, method = p_adj_method)
            n <- sum((d$coverage - d$alt_depth) >= alt_count_threshold & pval <= p_threshold, na.rm = TRUE)
            data.frame(loc = xi, mut_cell_n = n)
        }) %>% rbindlist()
        res <- res_2[mut_cell_n >= min_cell]
    }

    mtmutObj$loc_pass <- res$loc
    mtmutObj$loc_filter <- list(
        min_cell = min_cell,
        model = model,
        p_threshold = p_threshold,
        alt_count_threshold = alt_count_threshold,
        p_adj_method = p_adj_method
    )
    mtmutObj
}


#' Export the mutation matrix
#'
#' The helper functions to export the mutation results for further analysis. The output format can be data.frame, data.table or matrix for p value, allele frequency or binary mutation status.
#'
#' @param mtmutObj The scMtioMut object.
#' @param percent_interp A numeric value, the overlapping percentage threshold for triggering interpolation. The default is 1, which means no interpolation.
#' @param n_interp A integer value, the minimum number of overlapped cells with mutation for triggering interpolation, the default is 3.
#' @param all_cell A boolean to indicate whether to include all cells or only cells with mutation. By default, only cells with mutation are included.
#' @param memoSort A boolean to indicate whether to sort the loci by mutation frequency. The default is TRUE, the advanced user will know when to set it to FALSE.
#' @param ... Other parameters passed to \code{\link{export_dt}} or \code{\link{export_df}}.
#' @return data.frame, data.table or matrix or p
#' @export
#' @examples
#' ## Use the example data
#' f <- system.file("extdata", "mini_dataset.tsv.gz", package = "scMitoMut")
#'
#' ## Create a temporary h5 file
#' ## In real case, we keep the h5 in project folder for future use
#' f_h5_tmp <- tempfile(fileext = ".h5")
#'
#' ## Load the data with parse_table function
#' f_h5 <- parse_table(f, sep = "\t", h5_file = f_h5_tmp)
#' ## open the h5 file and create a mtmutObj object
#' x <- open_h5_file(f_h5)
#' run_model_fit(x)
#' x <- filter_loc(x, min_cell = 5, model = "bb", p_threshold = 0.05, p_adj_method = "fdr")
#' x
#' export_df(x)
#' export_pval(x)
#' export_af(x)
#' export_binary(x)
export_dt <- function(mtmutObj, percent_interp = 1, n_interp = 3, all_cell = FALSE) {

    if (!is(mtmutObj, "mtmutObj")) {
        stop("mtmutObj should be a mtmutObj object")
    }

    if (percent_interp < 0 | percent_interp > 1) {
        stop("percent_interp should be in [0, 1]")
    }

    if (n_interp < 0) {
        stop("n_interp should be >= 0")
    }

    if (all_cell != TRUE & all_cell != FALSE) {
        stop("all_cell should be TRUE or FALSE")
    }

    loc_list <- mtmutObj$loc_pass

    ## Get the parameters from mtmutObj
    min_cell <- mtmutObj$loc_filter$min_cell
    model <- mtmutObj$loc_filter$model
    p_threshold <- mtmutObj$loc_filter$p_threshold
    alt_count_threshold <- mtmutObj$loc_filter$alt_count_threshold
    p_adj_method <- mtmutObj$loc_filter$p_adj_method

    res <- parallel::mclapply(loc_list, function(loc_i) {
        d <- read_locus(mtmutObj, loc_i)
        d$pval <- get_pval(mtmutObj, loc_i, model = model, method = p_adj_method)
        d$loc <- loc_i
        d
    }) %>% rbindlist()

    res[, af := ((fwd_depth + rev_depth) / coverage)]
    res[, alt_count := (coverage - alt_depth)]

    # if (!all_cell) {
    # cell_list = res[, .(n = sum(pval < p_threshold, na.rm=TRUE)), by = cell_barcode][n > 0, cell_barcode]
    # loc_list = res[, .(n = sum(pval < p_threshold, na.rm=TRUE)), by = loc][n >= min_cell, loc]
    # res <- res[cell_barcode %in% cell_list & loc %in% loc_list]
    # }


    ## matrix pval
    d <- data.table::dcast(res, loc ~ cell_barcode, value.var = "pval")
    m_p <- data.matrix(d[, -1])
    rownames(m_p) <- d[[1]]

    ## matrix af
    d <- data.table::dcast(res, loc ~ cell_barcode, value.var = "af")
    m_a <- data.matrix(d[, -1])
    rownames(m_a) <- d[[1]]

    ## matrix alt_count
    d <- data.table::dcast(res, loc ~ cell_barcode, value.var = "alt_count")
    m_c <- data.matrix(d[, -1])
    rownames(m_c) <- d[[1]]

    ## matrix binary
    m_b <- m_p < p_threshold & m_c >= alt_count_threshold

    ## interpolate the mutation has higher frequency with the mutation has lower frequency
    if (percent_interp < 1) {
        mut_count <- rowSums(m_b, na.rm = TRUE)
        ix <- sort(mut_count, decreasing = FALSE, index.return = TRUE)$ix

        ## For the lowest framqent mutations
        for (i in seq_len(length(ix) - 1)) {
            mut_i <- ix[i]

            ## test the overlap with other mutations
            for (j in (i + 1):(length(ix))) {
                mut_j <- ix[j]
                mut_v_maj <- m_b[mut_j, ]
                mut_v_min <- m_b[mut_i, ]
                mut_v_maj[is.na(mut_v_maj)] <- FALSE
                mut_v_min[is.na(mut_v_min)] <- FALSE

                # p <- fisher.test(mut_v_maj, mut_v_min)$p.value
                tab <- table(mut_v_maj, mut_v_min)
                p <- (tab[4] / sum(tab[, 2]))

                ## if the p value is small, then the two mutations are not independent
                ## interpolate the mutation has higher frequency with the mutation has lower frequency
                if (p >= percent_interp & tab[4] >= n_interp) {
                    m_b[mut_j, mut_v_min] <- TRUE
                } else {
                    m_b[mut_j, mut_v_min] <- FALSE
                }
            }
        }
    }
    m_b_dt <- m_b %>%
        data.table(keep.rownames = TRUE) %>%
        data.table::melt(id.var = "rn", value.name = "mut_status")
    res <- merge(res, m_b_dt, by.x = c("loc", "cell_barcode"), by.y = c("rn", "variable"), all = TRUE)

    if (!all_cell) {
        cell_list <- res[, .(n = sum(mut_status, na.rm = TRUE)), by = cell_barcode][n > 0, cell_barcode]
        loc_list <- res[, .(n = sum(mut_status, na.rm = TRUE)), by = loc][n >= min_cell, loc]
        res <- res[cell_barcode %in% cell_list & loc %in% loc_list]
    }

    return(res)
}

#' @export
#' @rdname export_dt
export_df <- function(mtmutObj, ...) {
    data.frame(export_dt(mtmutObj, ...))
}

#' @export
#' @rdname export_dt
export_pval <- function(mtmutObj, memoSort = TRUE, ...) {
    res <- export_dt(mtmutObj, ...)
    d <- data.table::dcast(res, loc ~ cell_barcode, value.var = "pval")
    m <- data.matrix(d[, -1])
    rownames(m) <- d[[1]]

    if (memoSort) {
        d <- data.table::dcast(res, loc ~ cell_barcode, value.var = "mut_status")
        m_b <- data.matrix(d[, -1])
        rownames(m_b) <- d[[1]]

        m <- memoSort(m_b, m)
    }

    m
}

#' @export
#' @rdname export_dt
export_binary <- function(mtmutObj, memoSort = TRUE, ...) {
    res <- export_dt(mtmutObj, ...)
    d <- data.table::dcast(res, loc ~ cell_barcode, value.var = "mut_status")
    m_b <- data.matrix(d[, -1])
    rownames(m_b) <- d[[1]]

    if (memoSort) {
        m_b <- memoSort(m_b)
    }

    m_b
}

#' @export
#' @rdname export_dt
export_af <- function(mtmutObj, memoSort = TRUE, ...) {
    res <- export_dt(mtmutObj, ...)
    d <- data.table::dcast(res, loc ~ cell_barcode, value.var = "af")
    m <- data.matrix(d[, -1])
    rownames(m) <- d[[1]]

    if (memoSort) {
        d <- data.table::dcast(res, loc ~ cell_barcode, value.var = "mut_status")
        m_b <- data.matrix(d[, -1])
        rownames(m_b) <- d[[1]]

        m <- memoSort(m_b, m)
    }

    m
}
