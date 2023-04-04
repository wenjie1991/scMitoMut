#######################################################################
#                          Private functions                          #
#######################################################################
## TODO: use the code elsewhere
# if (F) {
#     par(mfrow = c(2, 2))
#     ## alt / depth
#     plot(N, y / (N + 0.001), log = "x", xlab = "Seq Depth", ylab = "p", main = paste0(loc_i, "_", maj_base))
#     points(N[p_adj < p_threshold], (y / (N + 0.001))[p_adj < p_threshold], col = "red", pch = 19)
#     ## fwd alter / depth
#     plot(N, y_fwd / (N + 0.001), log = "x", xlab = "Seq Depth", ylab = "p", main = paste0(loc_i, "_", maj_base))
#     points(N[p_adj < p_threshold], (y_fwd / (N + 0.001))[p_adj < p_threshold], col = "red", pch = 19)
#     ## rev alter / depth
#     plot(N, y_rev / (N + 0.001), log = "x", xlab = "Seq Depth", ylab = "p", main = paste0(loc_i, "_", maj_base))
#     points(N[p_adj < p_threshold], (y_rev / (N + 0.001))[p_adj < p_threshold], col = "red", pch = 19)
#     ## fwd / rev
#     plot(N, y_fwd / (y_rev + 0.001), log = "x", xlab = "Seq Depth", ylab = "p", main = paste0(loc_i, "_", maj_base))
#     points(N[p_adj < p_threshold], (y_fwd / (y_rev + 0.001))[p_adj < p_threshold], col = "red", pch = 19)
# }

## TODO: use it in the future
# if (F) {
#     plot(N, y / N, log = 'x', xlab = "Seq Depth", ylab = "p", main = paste0(loc_i, "_", alt_base))
#     points(N[i], (y / N)[i], col = 'red', pch = 19)
# }

#' Draw depth ~ AF scatter plot
#'
#'
#' @export
plot_locus <- function(
    d_select_maj_base, p, p_threshold = 0.05, loc = NA, maj_base= NA) {
    ## alt / depth
    N = d_select_maj_base$coverage
    y = d_select_maj_base[, alt_depth]
    d_plot = data.table(
        depth = N,
        af = y / (N + 0.001),
        highlight = p < p_threshold
    )
    ggplot(d_plot, aes(x = depth, y = af, color = highlight)) + 
        geom_point() + 
        scale_color_manual(values = c("black", "red")) + 
        scale_x_log10() + 
        labs(x = "Seq Depth", y = "p", title = paste0(loc))#, "_", maj_base))
    # plot(N, y / (N + 0.001), log = "x", xlab = "Seq Depth", ylab = "p", 
        # main = paste0(loc_i, "_", maj_base))
    # points(N[p< p_threshold], (y / (N + 0.001))[p < p_threshold], 
        # col = "red", pch = 19)
}


#' @export
plot_af_coverage = function(mtmutObj, loc, p_threshold = 0.05, p_adj_method = "fdr") {
    d = read_locus(mtmutObj, loc)
    plot_locus(d, p = get_pval(mtmutObj, loc, p_adj_method), p_threshold = p_threshold, loc = loc)
}

#' QC plot: 2D scatter plot for coverage ~ AF and UMAP
#'
#' @param mtmutObj an object of class "mtmutObj"
#' @param loc a string of genome location, e.g. "mt1000"
#' @param p_threshold a numeric value of p-value threshold
#' @param p_adj_method a string of p-value adjustment method
#' @export
#' @examples
#' ##
plot_locus_profile = function(loc, mtmutObj, seuratObj, p_threshold = 0.05, p_adj_method = "fdr") {
    y = process_locus_bmbb(loc, mtmutObj)

    p_depth_af = plot_locus(y$d, p = p.adjust(y$model$beta_binom$pval, p_adj_method), p_threshold = p_threshold)
    cell_select_bb = y[[1]]$cell_barcode[p.adjust(y$model$beta_binom$pval, p_adj_method) < p_threshold]

    p_umap = Seurat::DimPlot(seuratObj, cells.highlight = cell_select_bb)

    egg::ggarrange(p_depth_af, p_umap, ncol = 2)
}

#' Heatmap plot
#' 
#' @param mtmutObj an object of class "mtmutObj"
#' @param loc_list a vector of mt genome location
#' @param cell_ann a vector of cell annotation
#' @param ann_colors a vector of cell annotation colors
#' @param type a string of plot type, "p" for p-value, "af" for allele frequency
#' @param ... other parameters for \code{\link{export_df}}
#'
#' @export
plot_heatmap = function(mtmutObj, loc_list, cell_ann = NULL, ann_colors = NULL, type = "p", ...) {

    if (type == "p") {
        ## heatmap of p value
        m = export_pval(mtmutObj, loc_list, memoSort = T, ...)

        p = pheatmap::pheatmap(m, 
            color = colorRampPalette((RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
            show_colnames = F, annotation_col = cell_ann, cluster_cols = F, 
            cluster_rows = F, annotation_colors = ann_colors)

    } else if (type == "af") {
        ## heatmap of af
        m = export_af(mtmutObj, loc_list, memoSort = T, ...)

        p = pheatmap::pheatmap(m, color = colorRampPalette((RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
            show_colnames = F, annotation_col = cell_ann, cluster_cols = F, 
            cluster_rows = F, annotation_colors = ann_colors)

    } else if (type == "binary") {
        ## heatmap of binary mutation
        m_b = export_binary(mtmutObj, loc_list, memoSort = T, ...)
        m_b %<>% as.data.frame() %>% data.matrix()

        p = pheatmap::pheatmap(m_b, show_colnames = F, annotation_col = cell_ann, 
            annotation_colors = ann_colors, cluster_cols = F, cluster_rows = F, legend = F)

    } else {
        stop("type should be either p, af or binary")
    }
    return(p)
}

