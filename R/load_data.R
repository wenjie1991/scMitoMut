## read data from mtGATK

read_mgatk = function(mgatk_output_dir, prefix) {
    base_reads_files = dir(mgatk_output_dir, str_glue("{prefix}.[ACTG].txt.gz"), full=T)
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


#' Load mgatk 
#'
#' @examples
#' 
#' y = load_mgatk(mgatk_output_dir, prefix = "lib_id", h5_file = "./mut.h5")
#'
#' ## cell list
#' y$cell_id
#' 
#' ## test h5 file read out
#' h5f = H5Fopen(y$h5_file)
#' h5f
#' h5g = h5f&"mut_table"
#' h5g$x1600
#' h5f$cell_id_list
#' h5closeAll()
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
    d_ply(merge_d, "loc", function(x) {
        d_name = paste0("mt", x$loc[1])
        x$loc = NULL
        h5write(x, h5g, d_name)
    })

    ## cell list
    cell_id = merge_d$cell_barcode %>% unique
    selected_cell = cell_id
    h5write(cell_id, h5f, "cell_id_list")
    h5write(selected_cell, h5f, "selected_cell")

    #     H5Fclose(h5f)
    #     H5Gclose(h5g)
    #     h5closeAll()

    h5_file = h5_file
}

## test read out
# h5f = H5Fopen(h5_file)
# h5f
# h5g = h5f&"mut_table"
# h5g$x1600
# h5f$cell_id_list
# h5closeAll()

#' Open h5 file
#'
#' @export
open_h5_file <- function(h5_file) {
    h5f <- H5Fopen(h5_file)
    h5g <- h5f & "mut_table"
    loc_list <- h5ls(h5g, recursive=F)$name 
    cell_id_list <- h5f$cell_id_list
    selected_cell <- h5f$selected_cell
    ## we call the output as mtmutObj
    mtmutObj = list(h5f = h5f, 
        mut_table = h5g,
        loc_list = loc_list, 
        loc_selected = loc_list,
        cell_id_list = cell_id_list, 
        selected_cell = selected_cell
    )
    class(mtmutObj) = "mtmutObj"
    mtmutObj
}

#' @export
format.mtmutObj <- function(x, ...) {
    cat("h5 file: mtmutObj\n")
}

#' @export
print.mtmutObj <- function(x, ...) {
    cat("h5 file: mtmutObj\n")
}


#' @export
is.mtmutObj <- function(x) inherits(x, "mtmutObj")


#' Select cell
#'
#' @export
subset_cell <- function(mtmutObj, cell_id_list) {
    h5delete(mtmutObj$h5f, "selected_cell")
    h5write(cell_id_list, mtmutObj$h5f, "selected_cell")
    mtmutObj$selected_cell <- cell_id_list
    mtmutObj
}

#' Select loci
#'
#' @export
subset_loc <- function(mtmutObj, loc_list) {
    mtmutObj$loc_selected <- loc_list
    mtmutObj
}

#' Get pval for single locus
#'
#' @export
get_pval = function(mtmutObj, loc, method = "fdr") {
    p.adjust((mtmutObj$h5f & "pval" & loc)[], method)
}

#' Call mutation
#'
#' @export
mut_filter = function(
    mtmutObj, min_cell = 1, p_threshold = 0.05, p_adj_method = "fdr") {
    loc_list = mtmutObj$loc_selected
    res = mclapply(loc_list, function(xi) {
        pval = get_pval(mtmutObj, xi, method = p_adj_method)
        data.frame(loc = xi, mut_cell_n = sum(pval < p_threshold, na.rm=T))
    }) %>% rbindlist
    res[mut_cell_n >= min_cell]
}



#' Read one locus
read_locus = function(mtmutObj, loc, maj_base = NULL) {

    d_sub <- data.table((mtmutObj$mut_table&loc)[])[cell_barcode %in% mtmutObj$selected_cell]

    ## maj_base is the base with max frequency
    if (is.null(maj_base)) {
        # base_f <- table(d_sub$alt)
        # maj_base <- names(base_f)[which.max(base_f)]
        maj_base <- d_sub[, 
            .(vaf = median((fwd_depth + rev_depth) / coverage)), 
            by = alt][which.max(vaf), alt]
    }

    d_coverage <- unique(d_sub[, .(cell_barcode, coverage)])

    d_select_maj_base <- d_sub[alt == maj_base, .(cell_barcode, alt_depth = fwd_depth + rev_depth, fwd_depth, rev_depth)]
    d_select_maj_base <- merge(d_select_maj_base, d_coverage, by = "cell_barcode", all = T)
    d_select_maj_base[is.na(d_select_maj_base)] <- 0

    ## TODO: do we need this?
    # d_select_maj_base = d_select_maj_base[alt_depth / coverage > 0.8]

    d_select_maj_base
}


# This function sorts the matrix for better visualization of mutual exclusivity across genes
# Adaptoed from: https://gist.github.com/armish/564a65ab874a770e2c26
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



#' Export the mutation matrix
#' 
#' @param mtmutObj The scMtioMut object
#' @param pos_list The list of positions
#' @param p_threshold The p value threshold
#' @param min_cell_n The minimum number of cells with mutation
#' @param p_adj_method The method to adjust p value
#' @return The mutation matrix
#' @export
#' @examples
#' # d_plot = export_p(x, pos_list)
#' # d_plot_vaf = export_vaf(x, pos_list)
#' # d_plot_bin = export_binary(x, pos_list, 0.05)
export_dt = function(mtmutObj, pos_list, p_threshold = 0.05, percent_interp = 0.2, n_interp = 2, min_cell_n = 0, p_adj_method = "fdr", all_cell = F) {
    res = mclapply(pos_list, function(pos_i) {
        d = read_locus(mtmutObj, pos_i)
        d$pval = get_pval(mtmutObj, pos_i, method = "fdr")
        d$pos = pos_i
        d
    }) %>% rbindlist
    
    res[, vaf := ((fwd_depth + rev_depth) / coverage)]

    # if (!all_cell) {
        # cell_list = res[, .(n = sum(pval < p_threshold, na.rm=T)), by = cell_barcode][n > 0, cell_barcode]
        # pos_list = res[, .(n = sum(pval < p_threshold, na.rm=T)), by = pos][n >= min_cell_n, pos]
        # res = res[cell_barcode %in% cell_list & pos %in% pos_list]
    # }

    ## get binary mutation
    d = dcast(res, pos ~ cell_barcode, value.var = "pval")
    m = data.matrix(d[, -1])
    rownames(m) =  d[[1]]
    m_b = m < p_threshold

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
    res = merge(res, m_b_dt, by.x = c('pos', "cell_barcode"), by.y = c('rn', "variable"), all = T)

    if (!all_cell) {
        cell_list = res[, .(n = sum(mut_status, na.rm=T)), by = cell_barcode][n > 0, cell_barcode]
        pos_list = res[, .(n = sum(mut_status, na.rm=T)), by = pos][n >= min_cell_n, pos]
        res = res[cell_barcode %in% cell_list & pos %in% pos_list]
    }

    return(res)
}

export_df = function(mtmutObj, pos_list, ...) {
    data.frame(export_dt(mtmutObj, pos_list, ...))
}

export_pval = function(mtmutObj, pos_list, memoSort = F, ...) {
    res = export_dt(mtmutObj, pos_list, ...)
    d = dcast(res, pos ~ cell_barcode, value.var = "pval")
    m = data.matrix(d[, -1])
    rownames(m) =  d[[1]]

    if (memoSort) {
        d = dcast(res, pos ~ cell_barcode, value.var = "mut_status")
        m_b = data.matrix(d[, -1])
        rownames(m_b) =  d[[1]]

        m = memoSort(m_b, m)
    }

    m
}

export_binary = function(mtmutObj, pos_list, memoSort = F, ...) {
    res = export_dt(mtmutObj, pos_list, ...)
    d = dcast(res, pos ~ cell_barcode, value.var = "mut_status")
    m_b = data.matrix(d[, -1])
    rownames(m_b) =  d[[1]]

    if (memoSort) {
        m_b = memoSort(m_b)
    }

    m_b
}

export_vaf = function(mtmutObj, pos_list, memoSort = F, ...) {
    res = export_dt(mtmutObj, pos_list, ...)
    d = dcast(res, pos ~ cell_barcode, value.var = "vaf")
    m = data.matrix(d[, -1])
    rownames(m) =  d[[1]]

    if (memoSort) {
        d = dcast(res, pos ~ cell_barcode, value.var = "mut_status")
        m_b = data.matrix(d[, -1])
        rownames(m_b) =  d[[1]]

        m = memoSort(m_b, m)
    }

    m
}


