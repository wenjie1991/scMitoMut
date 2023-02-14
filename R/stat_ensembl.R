#' @export
process_locus_ensembl <- function(loc, mtmutObj, maj_base = NULL, max_try = 1) {
    d_select_maj_base <- read_locus(mtmutObj, loc, maj_base)

    res_bb <- process_locus_bb(d_select_maj_base, max_try = max_try)
    res_bm <- process_locus_bm(d_select_maj_base, max_try = max_try)

    res <- list(
        locus = d_select_maj_base,
        model = list(
            beta_binom = res_bb,
            binom_mix = res_bm
        )
    )
}
