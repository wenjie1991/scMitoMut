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
fit_bb <- function(y, N, init_vaf, max_try = 1) {
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

#' Calculate p-value
#'
calc_pval <- function(y, N, fit) {

    pval <- NULL

    if (class(fit) != "try-error") {

        prob <- fit[1]
        theta <- fit[2]

        pval <- lapply(seq_along(y), function(i) {
            sum(emdbook::dbetabinom(0:(y[i]), prob, N[i], theta))
        }) %>% unlist()
    }

    pval
}

#' Fit one location
#'
#' @examples
#'
#' process_locus("mt1000", mtmutObj, maj_base = NULL, max_try = 1)
#'
#' ##
process_locus_bb <- function(d_select_maj_base, max_try = 1) {

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
    fit = fit_bb(y, N, mean_vaf, max_try)

    pval <- calc_pval(y, N, fit)

    #################################################
    ## Output
    model_par <- data.table(mean_vaf, prob = fit[1], theta = fit[2])

    ## do p value adjust later
    ## and join the p value to the table later
    # d_select_maj_base$p <- pval

    ## output the list
    list(
        pval = pval,
        parameters = model_par
    )
}


