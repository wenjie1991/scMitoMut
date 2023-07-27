#######################################################################
#                         Internal function                           #
#######################################################################
get_bm_pval = function(x, method = "none") {
    p.adjust(x$bm_pval, method)
}

#######################################################################
#                   End of internal function region                   #
#######################################################################



#' Fit tree models for one locus
#' 
#' His is a wrapper function for \code{\link{process_locus_bm}} and \code{\link{process_locus_bb}}. 
#' It will fit binomial mixture model, beta binomial model and calculate the VMR and consistency of fwd rev strand.
#'
#' @param mtmutObj a mtmutObj object.
#' @param loc string given the locus name (e.g. "mt1000").
#' @param maj_base string given the dominant allele (e.g. "A"), if NULL auto detect the dominant allele.
#' @param return_data logical whether to return the allele count data, if FALSE, the \code{data} in the return value will be NULL. The default is FALSE.
#' @param ... other parameters passed to \code{\link{process_locus_bm}} and \code{\link{process_locus_bb}}.
#' @return A list of three elements:
#' \item{data}{data.frame of the allele count data.}
#' \item{locus}{data.table of the VMR and consistency of fwd rev strand.}
#' \item{model}{list of the model fitting results.}
#' @example
#' 
#' @export
process_locus_bmbb <- function(mtmutObj, loc, maj_base = NULL, return_data = FALSE, ...) {
    d_select_maj_base <- read_locus(mtmutObj, loc, maj_base)

    ## fit binomial mixture model
    res_bm <- process_locus_bm(d_select_maj_base, ...)

    ## fit beta binomial model
    selected_maj_cell = d_select_maj_base[get_bm_pval(res_bm, 'fdr') >= 0.05]$cell_barcode
    # selected_maj_cell = d_select_maj_base[get_bm_pval(res_bm, 'fdr') >= 0.01]$cell_barcode
    res_bb <- process_locus_bb(d_select_maj_base, selected_maj_cell, ...)

    ## VMR and consistency of fwd rev strand
    res_summary <- process_locus_summary(d_select_maj_base)

    if (!return_data) {
        d_select_maj_base = NULL
    }

    res <- list(
        data = d_select_maj_base,
        locus = res_summary,
        model = list(
            beta_binom = res_bb,
            binom_mix = res_bm
        )
    )
}


#' Fit binomial mixture model for every candidate locus
#'
#' @param mtmutObj a mtmutObj object.
#' @param mc.cores integer number of cores to use.
#' @return NULL, the results are saved in the h5f file.
#' @details
#' This function will fit three models for every candidate locus:
#' \itemize{
#' \item{binomial mixture model}{}
#' \item{beta binomial model}{}
#' \item{binomial model}{}
#' }
#' The results are saved in the h5f file.
#' @example
#' @export
run_model_fit <- function(mtmutObj, mc.cores = getOption("mc.cores", 2L)) {
    ## get the list of loci
    loc_list <- mtmutObj$loc_selected

    ## create a h5g to keep p value
    if ("pval" %in% h5ls(mtmutObj$h5f, recursive = F)$name) {
        h5delete(mtmutObj$h5f, "pval")
    }
    h5g <- H5Gcreate(h5loc = mtmutObj$h5f, name = "pval")

    ## run the ensemble calling
    # pb <- progress::progress_bar$new(total = length(loc_list))
    res_l = parallel::mclapply(loc_list, function(xi) {
        print(xi)
        # pb$tick()
        res = process_locus_bmbb(mtmutObj, xi, return_data = FALSE)

        list(
            bi_pval = res$model$binom_mix$bi_pval,
            bm_pval = res$model$binom_mix$bm_pval,
            bb_pval = res$model$beta_binom$bb_pval,
            model_par_bb = res$model$beta_binom$parameters,
            model_par_bm = res$model$binom_mix$parameters
        )
    }, mc.cores = mc.cores)

    model_par_bb = lapply(res_l, function(x) {
        x$model_par_bb
    }) %>% rbindlist #%>% do.call(rbind, .) %>% data.frame
    model_par_bb = cbind(loc = loc_list, model_par_bb)

    model_par_bm = lapply(res_l, function(x) {
        x$model_par_bm
    }) %>% rbindlist #%>% do.call(rbind, .) %>% data.frame
    model_par_bm = cbind(loc = loc_list, model_par_bm)

    lapply(seq_along(loc_list), function(i) {
        loc_i = loc_list[i]
        res_i = list(
            bi_pval = res_l[[i]]$bi_pval,
            bm_pval = res_l[[i]]$bm_pval,
            bb_pval = res_l[[i]]$bb_pval
        )
        h5write(res_i, h5g, loc_i)
    })

    if ("model_par_bb" %in% h5ls(x$h5f, recursive = F)$name) {
        h5delete(mtmutObj$h5f, "model_par_bb")
    }

    if ("model_par_bm" %in% h5ls(x$h5f, recursive = F)$name) {
        h5delete(mtmutObj$h5f, "model_par_bm")
    }

    h5write(model_par_bb, mtmutObj$h5f, "model_par_bb")
    h5write(model_par_bm, mtmutObj$h5f, "model_par_bm")
    gc()
}



