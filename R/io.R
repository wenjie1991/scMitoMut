######################################################################
#                          Internal function                          #
#######################################################################

## Parse the output of mtGATK into a data.table
read_mgatk = function(mgatk_output_dir, prefix) {
    mgatk_output_dir = paste0(mgatk_output_dir)
    base_reads_files = dir(mgatk_output_dir, str_glue("{prefix}.[ACTG].txt.gz"), full=T)
    print(base_reads_files)
    names(base_reads_files) = str_split_fixed(base_reads_files %>% basename, "\\.", 4)[, 2]
    ref_file = dir(mgatk_output_dir, "ref", full=T)
    coverage_file = dir(mgatk_output_dir, "coverage", full=T)

    ## base-wise sequence depth
    base_depth_d = lapply(seq_along(base_reads_files), function(i) {
        f = base_reads_files[i]
        d = fread(f)
        names(d) = c("loc", "cell_barcode", "fwd_depth", "rev_depth")
        d$alt = names(base_reads_files)[i]
        d
    }) %>% rbindlist

    ## reference sequence
    ref_d = fread(ref_file)
    names(ref_d) = c("loc", "ref")

    ## base-wise coverage
    if (length(coverage_file) == 0) {
        coverage_d = base_depth_d[, .(coverage = sum(fwd_depth + rev_depth, na.rm=T)), by = c("loc", "cell_barcode")]
    } else {
        coverage_d = fread(coverage_file)
        names(coverage_d) = c("loc", "cell_barcode", "coverage")
    }

    ## merge base-wise depth, ref seq, coverage
    merge_d = merge(base_depth_d, coverage_d, by = c("loc", "cell_barcode"))
    merge_d = merge(merge_d, ref_d, by = "loc")
    # rm(base_depth_d)
    # rm(coverage_d)
    # rm(ref_d)
    # gc()
    merge_d
}

## This function sorts the matrix for better visualization of mutual exclusivity across genes
## Adaptoed from: https://gist.github.com/armish/564a65ab874a770e2c26
memoSort <- function(M, m = NULL) {
	geneOrder <- sort(rowSums(M, na.rm=T), decreasing=TRUE, index.return=TRUE)$ix;
	scoreCol <- function(x) {
		score <- 0;
		for(i in 1:length(x)) {
			if(x[i] & !is.na(x[i])) {
				score <- score + 2^(length(x)-i);
			}
		}
		return(score);
	}
	scores <- apply(M[geneOrder, ], 2, scoreCol);
	sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
    if (is.null(m)) {
        return(M[geneOrder, sampleOrder]);
    } else {
        return(m[geneOrder, sampleOrder]);
    }
}

## Read one locus
read_locus = function(mtmutObj, loc, maj_base = NULL) {

    d_sub <- data.table((mtmutObj$mut_table&loc)[])[cell_barcode %in% mtmutObj$cell_selected]

    ## maj_base is the base with max frequency
    if (is.null(maj_base)) {
        # base_f <- table(d_sub$alt)
        # maj_base <- names(base_f)[which.max(base_f)]
        maj_base <- d_sub[, 
            .(af = median((fwd_depth + rev_depth) / coverage)), 
            by = alt][which.max(af), alt]
    }

    d_coverage <- unique(d_sub[, .(cell_barcode, coverage)])

    d_select_maj_base <- d_sub[alt == maj_base, .(cell_barcode, alt_depth = fwd_depth + rev_depth, fwd_depth, rev_depth)]
    d_select_maj_base <- merge(d_select_maj_base, d_coverage, by = "cell_barcode", all = T)
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
#' 
#' @export
parse_table = function(file, h5_file = 'mut.h5', ...) {
    merge_d = fread(file, ...)

    ##############################
    ## save to h5 file
    if (file.exists(h5_file)) {
        warning("\nH5 file exists, remove it.\n")
        file.remove(h5_file)
    }

    H5Fcreate(h5_file)
    h5f = H5Fopen(h5_file)

    ## mut_table
    h5g = H5Gcreate(h5loc = h5f, name = "mut_table")
    plyr::d_ply(merge_d, "loc", function(x) {
        d_name = paste0("mt", x$loc[1])
        x$loc = NULL
        h5write(x, h5g, d_name)
    })

    ## cell list
    cell_list = merge_d$cell_barcode %>% unique
    h5write(cell_list, h5f, "cell_list")
    h5write(cell_list, h5f, "cell_selected")

    ## loc list
    loc_list = merge_d$loc %>% unique %>% paste0("mt", .)
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


#' Load mtGATK output save to H5 file
#'
#' @param dir a string of the mtGATK output \code{final} fold
#' @param prefix a string of the prefix of the mtGATK output directory
#' @param h5_file a string of the output h5 file directory
#' @return a string of the output h5 file directory
#'
#' @examples
#' 
#' ## Not run
#' # y = load_mgatk("./mgatk_out/final/", prefix = "prefix_name", h5_file = "./mut.h5")
#' # print(y)
#' # > [1] './mut.h5',
#' loc_pass = c(),
#'
#' ##
#' @export
parse_mgatk = function(dir, prefix, h5_file = "mut.h5") {

    ##############################
    ## Read in data
    merge_d = read_mgatk(mgatk_output_ = dir, prefix = prefix)

    ##############################
    ## save to h5 file
    if (file.exists(h5_file)) {
        warning("\nH5 file exists, remove it.\n")
        file.remove(h5_file)
    }

    H5Fcreate(h5_file)
    h5f = H5Fopen(h5_file)

    ## mut_table
    h5g = H5Gcreate(h5loc = h5f, name = "mut_table")
    plyr::d_ply(merge_d, "loc", function(x) {
        d_name = paste0("mt", x$loc[1])
        x$loc = NULL
        h5write(x, h5g, d_name)
    })

    ## cell list
    cell_list = merge_d$cell_barcode %>% unique
    h5write(cell_list, h5f, "cell_list")
    h5write(cell_list, h5f, "cell_selected")

    ## loc list
    loc_list = merge_d$loc %>% unique %>% paste0("mt", .)
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
#' @param h5_file a string of the h5 file directory
#' @return a mtmutObj object
#'
#' @export
open_h5_file <- function(h5_file) {
    h5f <- H5Fopen(h5_file)
    h5g <- h5f & "mut_table"
    loc_list <- h5f$loc_list
    loc_selected = h5f$loc_selected
    cell_list <- h5f$cell_list
    cell_selected <- h5f$cell_selected
    ## we call the output as mtmutObj
    mtmutObj = list(h5f = h5f, 
        mut_table = h5g,
        loc_list = loc_list, 
        loc_selected = loc_selected,
        cell_list = cell_list, 
        cell_selected = cell_selected,
        loc_pass = c(),
        loc_filter = list(
            min_cell = 1,
            model = "bb",
            p_threshold = 0.05,
            p_adj_method = 'fdr',
            af_threshold = 1 
        )
    )
    class(mtmutObj) = "mtmutObj"
    mtmutObj
}

## TODO: make the H5 file S3 class

#' Mitochondiral mutation object
#'
#' A S3 class for handling mitochondrial mutation data
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{h5f}{H5 file handle}
#'   \item{mut_table}{allele count table H5 group handle}
#'   \item{loc_list}{list of available loci}
#'   \item{loc_selected}{selected loci}
#'   \item{cell_list}{list of available cell ids}
#'   \item{cell_selected}{selected cell ids}
#' }

#' @export
format.mtmutObj <- function(x) {
    cat("h5 file: mtmutObj\n")
}

#' @export
print.mtmutObj <- function(x) {
    cat("h5 file: mtmutObj\n")
}

#' @export
is.mtmutObj <- function(x) inherits(x, "mtmutObj")


#' Subset cell and loci
#' 
#' @param mtmutObj a mtmutObj object
#' @param cell_list a list of cell barcodes
#' @param loc_list a list of loci
#' @return a mtmutObj object with cell and loci selected
#'
#' @export
subset_cell <- function(mtmutObj, cell_list) {
    if ("cell_selected" %in% h5ls(mtmutObj$h5f, recursive = F)$name) {
        h5delete(mtmutObj$h5f, "cell_selected")
    }
    h5write(cell_list, mtmutObj$h5f, "cell_selected")
    mtmutObj$cell_selected <- cell_list
    mtmutObj
}

#' @rdname subset_cell
#' @export
subset_loc <- function(mtmutObj, loc_list) {
    if ("loc_selected" %in% h5ls(mtmutObj$h5f, recursive = F)$name) {
        h5delete(mtmutObj$h5f, "loc_selected")
    }
    h5write(loc_list, mtmutObj$h5f, "loc_selected")
    mtmutObj$loc_selected <- loc_list
    mtmutObj
}

#' Get p-value list for single locus
#'
#' @param mtmutObj a mtmutObj object
#' @param loc a string of the locus
#' @param method a string of the method for p-value adjustment
#'
#' @export
get_pval = function(mtmutObj, loc, model = "bi", method = "fdr") {
    if (model == "bb") {
        pval_item = "bb_pval"
    } else if (model == "bm") {
        pval_item = "bm_pval"
    } else if (model == "bi") {
        pval_item = "bi_pval"
    } else {
        stop("model should be bb or bi")
    }
    p.adjust((mtmutObj$h5f & "pval" & loc & pval_item)[], method)
}

#' Filter and visualize mutations
#'
#' @param mtmutObj a mtmutObj object
#' @param min_cell a integer of the minimum number of cells with mutation
#' @param model a string of the model for mutation calling, it can be "bb", "bm" or "bi"
#' @param p_threshold a numeric of the p-value threshold
#' @param p_adj_method a string of the method for p-value adjustment, 
#'   refer to \code{\link[p.adjust]{p.adjust}}
#' @param af_threshold a numeric of the majority allele frequency threshold,
#'   the major allele af < af_threshold will be considered as mutation.
#'   The default is 1, which means no filtering
#' @return a mtmutObj object with loci selected
#'
#' @export
mut_filter = function(
    mtmutObj, min_cell = NULL, model = NULL, p_threshold = NULL, p_adj_method = NULL, af_threshold = NULL) {
    ## Get the parameters from mtmutObj if they are NULL
    min_cell = ifelse(is.null(min_cell), mtmutObj$loc_filter$min_cell, min_cell)
    model = ifelse(is.null(model), mtmutObj$loc_filter$model, model)
    p_threshold = ifelse(is.null(p_threshold), mtmutObj$loc_filter$p_threshold, p_threshold)
    p_adj_method = ifelse(is.null(p_adj_method), mtmutObj$loc_filter$p_adj_method, p_adj_method)
    af_threshold = ifelse(is.null(af_threshold), mtmutObj$loc_filter$af_threshold, af_threshold)

    loc_list = mtmutObj$loc_selected
    res = parallel::mclapply(loc_list, function(xi) {
        pval = get_pval(mtmutObj, xi, model = model, method = p_adj_method)
        data.frame(loc = xi, mut_cell_n = sum(pval < p_threshold, na.rm=T))
    }) %>% rbindlist
    res = res[mut_cell_n >= min_cell]
    mtmutObj$loc_pass = res$loc 
    mtmutObj$loc_filter = list(
        min_cell = min_cell,
        model = model,
        p_threshold = p_threshold,
        p_adj_method = p_adj_method,
        af_threshold = af_threshold
    )
    mtmutObj
}


#' Export the mutation matrix
#' 
#' @param mtmutObj The scMtioMut object.
#' @param percent_interp A numeric value, the overlapping percentage threshold for triggering interpolation. The default is 1, which means no interpolation.
#' @param all_cell A boolean to indicate whether to include all cells or only cells with mutation.
#' @return data.frame, data.table or matrix
#' @export
#' @examples
#' # d_plot = export_p(x, loc_list)
#' # d_plot_af = export_af(x, loc_list)
#' # d_plot_bin = export_binary(x, loc_list, 0.05)
export_dt = function(mtmutObj, percent_interp = 1, n_interp = 3, all_cell = F) {

    loc_list = mtmutObj$loc_pass

    ## Get the parameters from mtmutObj
    min_cell = mtmutObj$loc_filter$min_cell
    model = mtmutObj$loc_filter$model
    p_threshold = mtmutObj$loc_filter$p_threshold
    p_adj_method = mtmutObj$loc_filter$p_adj_method
    af_threshold = mtmutObj$loc_filter$af_threshold

    res = parallel::mclapply(loc_list, function(loc_i) {
        d = read_locus(mtmutObj, loc_i)
        d$pval = get_pval(mtmutObj, loc_i, model = model, method = p_adj_method)
        d$loc = loc_i
        d
    }) %>% rbindlist
    
    res[, af := ((fwd_depth + rev_depth) / coverage)]

    # if (!all_cell) {
        # cell_list = res[, .(n = sum(pval < p_threshold, na.rm=T)), by = cell_barcode][n > 0, cell_barcode]
        # loc_list = res[, .(n = sum(pval < p_threshold, na.rm=T)), by = loc][n >= min_cell, loc]
        # res = res[cell_barcode %in% cell_list & loc %in% loc_list]
    # }


    ## matrix pval
    d = dcast(res, loc ~ cell_barcode, value.var = "pval")
    m_p = data.matrix(d[, -1])
    rownames(m_p) =  d[[1]]

    ## matrix af
    d = dcast(res, loc ~ cell_barcode, value.var = "af")
    m_a = data.matrix(d[, -1])
    rownames(m_a) =  d[[1]]

    ## matrix binary
    m_b = m_p < p_threshold & m_a < af_threshold

    ## interpolate the mutation has higher frequency with the mutation has lower frequency
    if (percent_interp < 1) {
        mut_count = rowSums(m_b, na.rm=T)
        ix = sort(mut_count, decreasing = F, index.return = T)$ix

        ## For the lowest framqent mutations
        for (i in 1:(length(ix) - 1)) {
            mut_i = ix[i]

            ## test the overlap with other mutations
            for (j in (i+1):(length(ix))) {
                mut_j = ix[j]
                mut_v_maj = m_b[mut_j, ]
                mut_v_min = m_b[mut_i, ]
                mut_v_maj[is.na(mut_v_maj)] = F
                mut_v_min[is.na(mut_v_min)] = F

                # p = fisher.test(mut_v_maj, mut_v_min)$p.value
                tab = table(mut_v_maj, mut_v_min)
                p = (tab[4] / sum(tab[, 2]))

                ## if the p value is small, then the two mutations are not independent
                ## interpolate the mutation has higher frequency with the mutation has lower frequency
                if (p >= percent_interp & tab[4] >= n_interp) {
                    m_b[mut_j, mut_v_min] = T
                } else {
                    m_b[mut_j, mut_v_min] = F
                }
            }
        }

    }
    m_b_dt = m_b %>% data.table(keep.rownames=T) %>% melt(id.var = 'rn', value.name = "mut_status") 
    res = merge(res, m_b_dt, by.x = c('loc', "cell_barcode"), by.y = c('rn', "variable"), all = T)

    if (!all_cell) {
        cell_list = res[, .(n = sum(mut_status, na.rm=T)), by = cell_barcode][n > 0, cell_barcode]
        loc_list = res[, .(n = sum(mut_status, na.rm=T)), by = loc][n >= min_cell, loc]
        res = res[cell_barcode %in% cell_list & loc %in% loc_list]
    }

    return(res)
}

#' @export
#' @rdname export_dt
export_df = function(mtmutObj, ...) {
    data.frame(export_dt(mtmutObj, ...))
}

#' @export
#' @rdname export_dt
export_pval = function(mtmutObj, memoSort = F, ...) {
    res = export_dt(mtmutObj, ...)
    d = dcast(res, loc ~ cell_barcode, value.var = "pval")
    m = data.matrix(d[, -1])
    rownames(m) =  d[[1]]

    if (memoSort) {
        d = dcast(res, loc ~ cell_barcode, value.var = "mut_status")
        m_b = data.matrix(d[, -1])
        rownames(m_b) =  d[[1]]

        m = memoSort(m_b, m)
    }

    m
}

#' @export
#' @rdname export_dt
export_binary = function(mtmutObj, memoSort = F, ...) {
    res = export_dt(mtmutObj, ...)
    d = dcast(res, loc ~ cell_barcode, value.var = "mut_status")
    m_b = data.matrix(d[, -1])
    rownames(m_b) =  d[[1]]

    if (memoSort) {
        m_b = memoSort(m_b)
    }

    m_b
}

#' @export
#' @rdname export_dt
export_af = function(mtmutObj, memoSort = F, ...) {
    res = export_dt(mtmutObj, ...)
    d = dcast(res, loc ~ cell_barcode, value.var = "af")
    m = data.matrix(d[, -1])
    rownames(m) =  d[[1]]

    if (memoSort) {
        d = dcast(res, loc ~ cell_barcode, value.var = "mut_status")
        m_b = data.matrix(d[, -1])
        rownames(m_b) =  d[[1]]

        m = memoSort(m_b, m)
    }

    m
}


