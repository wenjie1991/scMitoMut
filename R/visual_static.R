#######################################################################
#                         Internal functions                          #
#######################################################################
## TODO: use the code elsewhere
# if (FALSE) {
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
# if (FALSE) {
#     plot(N, y / N, log = 'x', xlab = "Seq Depth", ylab = "p", main = paste0(loc_i, "_", alt_base))
#     points(N[i], (y / N)[i], col = 'red', pch = 19)
# }

# Draw depth ~ AF scatter plot
plot_locus <- function(
    d_select_maj_base, p, p_threshold = 0.05, loc = NA, maj_base = NA) {
  ## alt / depth
  N <- d_select_maj_base$coverage
  y <- d_select_maj_base[, alt_depth]
  d_plot <- data.table(
    depth = N,
    af = y / (N + 0.001),
    highlight = p < p_threshold
  )
  ggplot(d_plot, aes(x = depth, y = af, color = highlight)) +
    geom_point() +
    scale_color_manual(values = c("black", "red")) +
    scale_x_log10() +
    labs(x = "Depth", y = "Allele Frequency", color = "Mutant") +
    theme_bw()
  # , title = paste0(loc))#, "_", maj_base))
  # plot(N, y / (N + 0.001), log = "x", xlab = "Seq Depth", ylab = "p",
  # main = paste0(loc_i, "_", maj_base))
  # points(N[p< p_threshold], (y / (N + 0.001))[p < p_threshold],
  # col = "red", pch = 19)
}

# QC plot: 2D scatter plot for coverage ~ AF and UMAP
#
# @param loc a string of genome location, e.g. "chrM.200".
# @param seuratObj an object of class "Seurat".
# @param model a string of model name, one of "bb", "bm", "bi".
# @param p_threshold a numeric value of p-value threshold.
# @param p_adj_method a string of p-value adjustment method.
plot_locus_profile <- function(mtmutObj, loc, seuratObj, model = NULL, p_threshold = NULL, p_adj_method = NULL) {
  ## get parameters
  model <- ifelse(is.null(model), mtmutObj$loc_filter$model, model)
  p_threshold <- ifelse(is.null(p_threshold), mtmutObj$loc_filter$p_threshold, p_threshold)
  p_adj_method <- ifelse(is.null(p_adj_method), mtmutObj$loc_filter$p_adj_method, p_adj_method)

  y <- process_locus_bmbb(mtmutObj, loc, return_data = TRUE)
  if (model == "bb") {
    p_depth_af <- plot_locus(y$data, p = p.adjust(y$model$beta_binom$bb_pval, p_adj_method), p_threshold = p_threshold)
    cell_select_bb <- y$data$cell_barcode[p.adjust(y$model$beta_binom$bb_pval, p_adj_method) < p_threshold]
  } else if (model == "bm") {
    p_depth_af <- plot_locus(y$data, p = p.adjust(y$model$binom_mix$bm_pval, p_adj_method), p_threshold = p_threshold)
    cell_select_bb <- y$data$cell_barcode[p.adjust(y$model$binom_mix$bm_pval, p_adj_method) < p_threshold]
  } else if (model == "bi") {
    p_depth_af <- plot_locus(y$data, p = p.adjust(y$model$binom_mix$bi_pval, p_adj_method), p_threshold = p_threshold)
    cell_select_bb <- y$data$cell_barcode[p.adjust(y$model$binom_mix$bi_pval, p_adj_method) < p_threshold]
  } else {
    stop("model should be one of bb, bm, bi")
  }
  p_umap <- Seurat::DimPlot(seuratObj, cells.highlight = cell_select_bb) + theme_bw() + NoLegend()
  egg::ggarrange(p_depth_af, p_umap, ncol = 2)
}



#######################################################################
#                   End of internal function region                   #
#######################################################################


#' QC plot: 2D scatter plot for coverage ~ AF
#'
#' @param mtmutObj an object of class "mtmutObj".
#' @param loc a string of genome location, e.g. "chrM.200".
#' @param model a string of model name, one of "bb", "bm", "bi".
#' @param p_threshold a numeric value of p-value threshold.
#' @param p_adj_method a string of p-value adjustment method.
#' @examples
#' ## Use the example data
#' f <- system.file("extdata", "mini_dataset.tsv.gz", package = "scMitoMut")
#'
#' ## Create a temporary h5 file
#' ## In real case, we keep the h5 in project folder for future use
#' f_h5_tmp <- tempfile(fileext = ".h5")
#'
#' ## Load the data with parse_table function
#' f_h5 <- parse_table(f, sep = "\t", h5_file = f_h5_tmp)
#'
#' # open the h5f file
#' x <- open_h5_file(f_h5)
#'
#' # run the model fit
#' run_model_fit(x)
#' x
#'
#' # Filter the loci based on the model fit results
#' x <- filter_loc(x, min_cell = 5, model = "bb", p_threshold = 0.05, p_adj_method = "fdr")
#'
#' # plot the locus profile for chrM.200
#' plot_af_coverage(x, "chrM.204")
#' @export
plot_af_coverage <- function(mtmutObj, loc, model = NULL, p_threshold = NULL, p_adj_method = NULL) {
  ## get parameters
  model <- ifelse(is.null(model), mtmutObj$loc_filter$model, model)
  p_threshold <- ifelse(is.null(p_threshold), mtmutObj$loc_filter$p_threshold, p_threshold)
  p_adj_method <- ifelse(is.null(p_adj_method), mtmutObj$loc_filter$p_adj_method, p_adj_method)

  d <- read_locus(mtmutObj, loc)
  plot_locus(d, p = get_pval(mtmutObj, loc, model, p_adj_method), p_threshold = p_threshold, loc = loc)
}

#' Heatmap plot
#'
#' @param mtmutObj an object of class "mtmutObj".
#' @param type a string of plot type, "p" for p-value, "af" for allele frequency.
#' @param cell_ann a data.frame of cell annotation, with rownames as cell barcodes, please refer to \code{\link[pheatmap]{pheatmap}} for details.
#' @param ann_colors a list of colors for cell annotation with cell annotation as names, please refer to \code{\link[pheatmap]{pheatmap}} for details.
#' @param ... other parameters for \code{\link[scMitoMut]{export_df}} and \code{\link[pheatmap]{pheatmap}}.
#' @examples
#' # load the data
#' ## Use the example data
#' f <- system.file("extdata", "mini_dataset.tsv.gz", package = "scMitoMut")
#'
#' ## Create a temporary h5 file
#' ## In real case, we keep the h5 in project folder for future use
#' f_h5_tmp <- tempfile(fileext = ".h5")
#'
#' ## Load the data with parse_table function
#' f_h5 <- parse_table(f, sep = "\t", h5_file = f_h5_tmp)
#'
#' # open the h5f file
#' x <- open_h5_file(f_h5)
#' # run the model fit
#' run_model_fit(x)
#' x
#' # Filter the loci based on the model fit results
#' x <- filter_loc(x, min_cell = 5, model = "bb", p_threshold = 0.05, p_adj_method = "fdr")
#'
#' # set the cell annotation
#' f <- system.file("extdata", "mini_dataset_cell_ann.csv", package = "scMitoMut")
#' cell_ann <- read.csv(f, row.names = 1)
#' # Prepare the color for cell annotation
#' colors <- c(
#'   "Cancer Epi" = "#f28482",
#'   Blood = "#f6bd60"
#' )
#' ann_colors <- list("SeuratCellTypes" = colors)
#'
#' # plot the heatmap for p-value
#' plot_heatmap(x, type = "p", cell_ann = cell_ann, ann_colors = ann_colors, percent_interp = 0.2)
#' # plot the heatmap for allele frequency
#' plot_heatmap(x, type = "af", cell_ann = cell_ann, ann_colors = ann_colors, percent_interp = 0.2)
#' # plot the heatmap for binary mutation
#' plot_heatmap(x, type = "binary", cell_ann = cell_ann, ann_colors = ann_colors, percent_interp = 0.2)
#' @export
plot_heatmap <- function(mtmutObj, type = "p", cell_ann = NULL, ann_colors = NULL, ...) {
  if (type == "p") {
    ## heatmap of p value
    m <- export_pval(mtmutObj, memoSort = TRUE, ...)

    p <- pheatmap::pheatmap(m,
      color = rev(colorRampPalette((RColorBrewer::brewer.pal(n = 7, name = "GnBu")))(100)),
      show_colnames = FALSE, annotation_col = cell_ann, cluster_cols = FALSE,
      cluster_rows = FALSE, annotation_colors = ann_colors, na_col = "white"
    )
  } else if (type == "af") {
    ## heatmap of af
    m <- export_af(mtmutObj, memoSort = TRUE, ...)

    p <- pheatmap::pheatmap(m,
      color = rev(colorRampPalette((RColorBrewer::brewer.pal(n = 7, name = "GnBu")))(100)),
      show_colnames = FALSE, annotation_col = cell_ann, cluster_cols = FALSE,
      cluster_rows = FALSE, annotation_colors = ann_colors, na_col = "white"
    )
  } else if (type == "binary") {
    ## heatmap of binary mutation
    m_b <- export_binary(mtmutObj, memoSort = TRUE, ...)
    m_b %<>% as.data.frame() %>% data.matrix()

    p <- pheatmap::pheatmap(m_b,
      color = colorRampPalette((RColorBrewer::brewer.pal(n = 7, name = "GnBu")))(100),
      show_colnames = FALSE, annotation_col = cell_ann,
      annotation_colors = ann_colors, cluster_cols = FALSE, cluster_rows = FALSE, legend = FALSE,
      na_col = "white", ...
    )
  } else {
    stop("type should be either p, af or binary")
  }
  return(p)
}
