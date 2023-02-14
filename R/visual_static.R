## TODO: use the code elsewhere
if (F) {
    par(mfrow = c(2, 2))
    ## alt / depth
    plot(N, y / (N + 0.001), log = "x", xlab = "Seq Depth", ylab = "p", main = paste0(loc_i, "_", maj_base))
    points(N[p_adj < p_threshold], (y / (N + 0.001))[p_adj < p_threshold], col = "red", pch = 19)
    ## fwd alter / depth
    plot(N, y_fwd / (N + 0.001), log = "x", xlab = "Seq Depth", ylab = "p", main = paste0(loc_i, "_", maj_base))
    points(N[p_adj < p_threshold], (y_fwd / (N + 0.001))[p_adj < p_threshold], col = "red", pch = 19)
    ## rev alter / depth
    plot(N, y_rev / (N + 0.001), log = "x", xlab = "Seq Depth", ylab = "p", main = paste0(loc_i, "_", maj_base))
    points(N[p_adj < p_threshold], (y_rev / (N + 0.001))[p_adj < p_threshold], col = "red", pch = 19)
    ## fwd / rev
    plot(N, y_fwd / (y_rev + 0.001), log = "x", xlab = "Seq Depth", ylab = "p", main = paste0(loc_i, "_", maj_base))
    points(N[p_adj < p_threshold], (y_fwd / (y_rev + 0.001))[p_adj < p_threshold], col = "red", pch = 19)
}

## TODO: use it in the future
if (F) {
    plot(N, y / N, log = 'x', xlab = "Seq Depth", ylab = "p", main = paste0(loc_i, "_", alt_base))
    points(N[i], (y / N)[i], col = 'red', pch = 19)
}

#' @export
plot_locus <- function(
    d_select_maj_base, p_adj = 0.05, p_threshold, loc_i = NA, maj_base= NA) {
    ## alt / depth
    N = d_select_maj_base$coverage
    y = d_select_maj_base[, alt_depth]
    plot(N, y / (N + 0.001), log = "x", xlab = "Seq Depth", ylab = "p", 
        main = paste0(loc_i, "_", maj_base))
    points(N[p_adj < p_threshold], (y / (N + 0.001))[p_adj < p_threshold], 
        col = "red", pch = 19)
}


