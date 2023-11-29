#######################################################################
#                         Internal function                           #
#######################################################################

# Fit betaBinomial distribution
fit_bb_cpp <- function(y, N, max_iter = 100, tol = 1e-6) {
    mle_bb(y, N, max_iter, tol)
}

calc_pval_cpp <- function(y, N, fit) {
    a <- fit[[1]]
    b <- fit[[2]]
    pbetabinom(y, N, a, b)
}

# Fit beta-binomial distribution for one locus
#
# @param d_select_maj_base data.frame of one locus
# @param selected_maj_cell vector of selected cell
# @param max_theta maximum theta
# @param max_iter maximum iteration
# @param tol tolerance of log likelihood to stop iteration
# @return list of p-value and model parameters
process_locus_bb <- function(
    d_select_maj_base,
    selected_maj_cell = NULL,
    max_theta = 1e6,
    max_iter = 100,
    tol = 1e-3) {
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
    ## average af
    mean_af <- sum(y) / sum(N)

    if (!is.null(selected_maj_cell)) {
        ## vector of selected cell
        is_selected_cell <- d_select_maj_base$cell_barcode %in% selected_maj_cell

        #################################################
        ## Model fitting
        # fit = vglm(cbind(y, N - y) ~ 1, family = betabinomialff, trace = TRUE, weights = rep(1, length(N)))
        # fit = fit_bb_r(y, N, mean_af, max_try)
        fit <- fit_bb_cpp(y[is_selected_cell], N[is_selected_cell], max_iter = max_iter, tol = tol)
    } else {
        #################################################
        ## Model fitting
        # fit = vglm(cbind(y, N - y) ~ 1, family = betabinomialff, trace = TRUE, weights = rep(1, length(N)))
        # fit = fit_bb_r(y, N, mean_af, max_try)
        fit <- fit_bb_cpp(y, N, tol = tol, max_iter = max_iter)
    }

    ## correct fit to avoid overflow
    if (!is.na(fit[[1]]) & (fit[[1]] + fit[[2]] > max_theta)) {
        scale_factor <- max_theta / (fit[[1]] + fit[[2]])
        fit[[1]] <- fit[[1]] * scale_factor
        fit[[2]] <- fit[[2]] * scale_factor
    }

    ## the prob and theta
    prob <- fit[[1]] / (fit[[1]] + fit[[2]])
    theta <- fit[[1]] + fit[[2]]

    ## if the bb model prob is too small, set pval to 1
    if (is.na(prob)) {
        pval <- rep(NA, length(y))
        # } else if (1 - prob < 1e-6) {
        # pval = rep(1, length(y))
    } else {
        pval <- calc_pval_cpp(y, N, fit)
    }

    #################################################
    ## Output
    model_par <- data.table(mean_af, prob = prob, theta = theta)

    ## if the af > expected prob, set pval to 1
    pval[y / N > prob] <- 1

    ## output the list
    list(
        bb_pval = pval,
        parameters = model_par
    )
}

#######################################################################
#                   End of internal function region                   #
#######################################################################
