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
        coverage_d = base_depth_d[, .(coverage = sum(fwd_depth + rev_depth)), by = c("loc", "cell_barcode")]
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
    h5g = H5Gcreate(h5loc = h5f, name = "mut_table")
    # h5delete(h5f, "mut_table")

    d_ply(merge_d, "loc", function(x) {
        d_name = paste0("mt", x$loc[1])
        #         print(d_name)
        x$loc = NULL
        h5write(x, h5g, d_name)
    })

    ## cell list
    cell_id = merge_d$cell_barcode %>% unique
    h5write(cell_id, h5f, "cell_id_list")

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
    mut_list <- h5ls(h5g)$name 
    cell_id_list <- h5f$cell_id_list
    ## we call the output as mtmutObj
    list(h5f = h5f, 
        mut_table = h5g, 
        mut_list = mut_list, 
        cell_id_list = cell_id_list, 
        selected_loc = cell_id_list
    )
}

#' Select cell
#'
#' @export
subset_cell <- function(mtmutObj, cell_id_list) {
    mtmutObj$selected_loc <- cell_id_list
    mtmutObj
}


#' Read one locus
read_locus = function(mtmutObj, loc, maj_base = NULL) {

    d_sub <- data.table((mtmutObj$mut_table&loc)[])[cell_barcode %in% mtmutObj$selected_loc]

    ## maj_base is the base with max frequency
    if (is.null(maj_base)) {
        base_f <- table(d_sub$alt)
        maj_base <- names(base_f)[which.max(base_f)]
    }

    d_coverage <- unique(d_sub[, .(cell_barcode, coverage)])

    d_select_maj_base <- d_sub[alt == maj_base, .(cell_barcode, alt_depth = fwd_depth + rev_depth, fwd_depth, rev_depth)]
    d_select_maj_base <- merge(d_select_maj_base, d_coverage, by = "cell_barcode", all = T)
    d_select_maj_base[is.na(d_select_maj_base)] <- 0

    ## TODO: do we need this?
    # d_select_maj_base = d_select_maj_base[alt_depth / coverage > 0.8]

    d_select_maj_base
}
