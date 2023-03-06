
get_bm_pval = function(x, method = "none") {
    p.adjust(x$pval, method)
}


#' @export
process_locus_ensemble <- function(loc, mtmutObj, maj_base = NULL, max_try = 1) {
    d_select_maj_base <- read_locus(mtmutObj, loc, maj_base)

    ## fit binomial mixture model
    res_bm <- process_locus_bm(d_select_maj_base, max_try = max_try)

    ## fit beta binomial model
    selected_maj_cell = d_select_maj_base[get_bm_pval(res_bm, 'fdr') >= 0.001]$cell_barcode
    res_bb <- process_locus_bb(d_select_maj_base, selected_maj_cell)

    ## VMR and consistency of fwd rev strand
    res_summary <- process_locus_summary(d_select_maj_base)

    res <- list(
        data = d_select_maj_base,
        locus = res_summary,
        model = list(
            beta_binom = res_bb,
            binom_mix = res_bm
        )
    )
}
