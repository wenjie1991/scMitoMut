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
    d_select_maj_base, p, p_threshold = 0.05, loc_i = NA, maj_base= NA) {
    ## alt / depth
    N = d_select_maj_base$coverage
    y = d_select_maj_base[, alt_depth]
    d_plot = data.table(
        depth = N,
        vaf = y / (N + 0.001),
        highlight = p < p_threshold
    )
    ggplot(d_plot, aes(x = depth, y = vaf, color = highlight)) + 
        geom_point() + 
        scale_color_manual(values = c("black", "red")) + 
        scale_x_log10() + 
        labs(x = "Seq Depth", y = "p", title = paste0(loc_i, "_", maj_base))
    # plot(N, y / (N + 0.001), log = "x", xlab = "Seq Depth", ylab = "p", 
        # main = paste0(loc_i, "_", maj_base))
    # points(N[p< p_threshold], (y / (N + 0.001))[p < p_threshold], 
        # col = "red", pch = 19)
}


#' @export
plot_vaf_coverage = function(mtmutObj, loc, p_threshold = 0.05, p_adj_method = "fdr") {
    d = read_locus(mtmutObj, loc)
    plot_locus(d, p = get_pval(mtmutObj, loc, p_adj_method), p_threshold = p_threshold)
}

#' QC plot: 2D scatter plot for coverage ~ VAF and UMAP
#'
#' @param mtmutObj an object of class "mtmutObj"
#' @param loc a string of mitochondrial position, e.g. "mt1000"
#' @param p_threshold a numeric value of p-value threshold
#' @param p_adj_method a string of p-value adjustment method
#' @export
#' @examples
#' try_plot("chr1_10000", mtmutObj, p = 0.05, p_adj_method = "holm")
#' ##
plot_locus_profile = function(loc, mtmutObj, seuratObj, p_threshold = 0.05, p_adj_method = "holm") {
    y = process_locus_ensemble(loc, mtmutObj)

    p_depth_vaf = plot_locus(y$d, p = p.adjust(y$model$beta_binom$pval, p_adj_method), p_threshold = p_threshold)
    cell_select_bb = y[[1]]$cell_barcode[p.adjust(y$model$beta_binom$pval, p_adj_method) < p_threshold]

    p_umap = Seurat::DimPlot(seuratObj, cells.highlight = cell_select_bb)

    egg::ggarrange(p_depth_vaf, p_umap, ncol = 2)
}
