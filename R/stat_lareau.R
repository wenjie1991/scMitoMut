
#' Calculate strand bias
calc_strand_concordance <- function(y_fwd, y_rev) {
    cor(y_fwd, y_rev)
}

#' Calculate AF variant
calc_vmr <- function(y, N) {
    af_v = y / N
    bulk_af = mean(af_v)

    var(af_v) / (bulk_af + 1e-12)
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
    ## averrage af
    mean_af <- sum(y) / sum(N)

    cor = calc_strand_concordance(y_fwd, y_rev)
    vmr = calc_vmr(y, N)

    res <- list(
        mean_af = mean_af,
        cor = cor,
        vmr = vmr
    )
}


