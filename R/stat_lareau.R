
#' Calculate strand bias
calc_strand_concordance <- function(y_fwd, y_rev) {
    cor(y_fwd, y_rev)
}

#' Calculate VAF variant
calc_vmr <- function(y, N) {
    vaf_v = y / N
    bulk_vaf = mean(vaf_v)

    var(vaf_v) / (bulk_vaf + 1e-12)
}

process_locus_summary <- function(d_select_maj_base) {

    #################################################
    ### Transform data
    ## coverage
    N <- d_select_maj_base$coverage
    ## majority base depth
    y <- d_select_maj_base[, alt_depth]
    ## majority base forward depth
    y_fwd <- d_select_maj_base[, fwd_depth]
    ## majority base reverse depth
    y_rev <- d_select_maj_base[, rev_depth]
    ## averrage vaf
    mean_vaf <- sum(y) / sum(N)

    cor = calc_strand_concordance(y_fwd, y_rev)
    vmr = calc_vmr(y, N)

    res <- list(
        mean_vaf = mean_vaf,
        cor = cor,
        vmr = vmr
    )
}


