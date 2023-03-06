# library(data.table)
# library(ggplot2)
# library(plyr)
# library(magrittr)
# library(readr)
# library(stringr)
# library(pheatmap)
# library(rhdf5)
# library(qvalue)

## mutation preprocess mgatk output to HDR5
library(data.table)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(rhdf5)

# No import
library(emdbook)
library(bbmle)



#' Fit betabinomial distribution
betaBinomEst <- function(counts, total, prob_start = 0.99, theta_start = 10) {
    mtmp <- function(prob, size, theta) {
        -sum(emdbook::dbetabinom(counts, prob, size, theta, log = TRUE))
    }
    m0 <- suppressWarnings(
        bbmle::mle2(mtmp, start = list(prob = prob_start, theta = theta_start), data = list(size = total))
    )
    # MLE of theta
    bbmle::coef(m0)
}

#' Fit betaBinomial distribution
#'
fit_bb_r <- function(y, N, init_vaf, max_try = 1) {
    fit <- NULL
    try_num <- 0
    while (try_num < max_try) {
        try_num <- try_num + 1
        tryCatch(
            {
                fit <- betaBinomEst(y, N, init_vaf)
                break
            },
            error = function(e) {
                print(e)
                print(paste0("Try ", try_num, " times."))
            }
        )
        init_vaf <- init_vaf - 0.001
    }
    fit
}

fit_bb_cpp <- function(y, N, max_iter = 100, tol = 1e-6) {
    mle_bb(y, N, max_iter, tol)
}


#' Calculate p-value
#'
calc_pval_r <- function(y, N, fit) {

    pval <- NULL

    if (class(fit) != "try-error") {

        a = fit[[1]]
        b = fit[[2]]

        prob <- a / (a + b)
        theta <- a + b

        pval <- lapply(seq_along(y), function(i) {
            sum(emdbook::dbetabinom(0:(y[i]), prob, N[i], theta))
        }) %>% unlist()
    }

    pval
}

calc_pval_cpp = function(y, N, fit) {
    a = fit[[1]]
    b = fit[[2]]
    pbetabinom(y, N, a, b)
}

#' Fit one location
#'
#' @examples
#'
#' process_locus("mt1000", mtmutObj, maj_base = NULL, max_try = 1)
#'
#' ##
process_locus_bb_bk <- function(d_select_maj_base, max_try = 1) {

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

    #################################################
    ## Model fitting
    # fit = vglm(cbind(y, N - y) ~ 1, family = betabinomialff, trace = T, weights = rep(1, length(N)))
	# fit = fit_bb_r(y, N, mean_vaf, max_try)
    fit = fit_bb_cpp(y, N)

    # pval <- calc_pval_r(y, N, fit)
    pval <- calc_pval_cpp(y, N, fit)

    #################################################
    ## Output
    prob = fit[[1]] / (fit[[1]] + fit[[2]])
    theta = fit[[1]] + fit[[2]]

    model_par <- data.table(mean_vaf, prob = prob, theta = theta)

    ## do p value adjust later
    ## and join the p value to the table later
    # d_select_maj_base$p <- pval

    ## output the list
    list(
        pval = pval,
        parameters = model_par
    )
}


## new fit bb function
process_locus_bb = function(d_select_maj_base, selected_maj_cell = NULL, max_theta = 1e4) {
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

    if (!is.null(selected_maj_cell)) {
        ## vector of selected cell
        is_selected_cell = d_select_maj_base$cell_barcode %in% selected_maj_cell

        #################################################
        ## Model fitting
        # fit = vglm(cbind(y, N - y) ~ 1, family = betabinomialff, trace = T, weights = rep(1, length(N)))
        # fit = fit_bb_r(y, N, mean_vaf, max_try)
        fit = fit_bb_cpp(y[is_selected_cell], N[is_selected_cell])
    } else {

        #################################################
        ## Model fitting
        # fit = vglm(cbind(y, N - y) ~ 1, family = betabinomialff, trace = T, weights = rep(1, length(N)))
        # fit = fit_bb_r(y, N, mean_vaf, max_try)
        fit = fit_bb_cpp(y, N)
    }

    ## correct fit to avoid overflow
    if (fit[[1]] + fit[[2]] > max_theta) {
        scale_factor = max_theta / (fit[[1]] + fit[[2]])
        fit[[1]] = fit[[1]] * scale_factor
        fit[[2]] = fit[[2]] * scale_factor
    }


    # pval <- calc_pval_r(y, N, fit)
    pval <- calc_pval_cpp(y, N, fit)

    #################################################
    ## Output
    prob = fit[[1]] / (fit[[1]] + fit[[2]])
    theta = fit[[1]] + fit[[2]]

    model_par <- data.table(mean_vaf, prob = prob, theta = theta)

    ## output the list
    list(
        pval = pval,
        parameters = model_par
    )
}

