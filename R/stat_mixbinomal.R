#######################################################################
#                         Internal functions                          #
#######################################################################
# fit_bm_r <- function(x, n, init_theta, init_lambda) {
#
#     m = data.matrix(data.frame(x, n-x))
#     em_out1 <- mixtools::multmixEM(m, k = 1)
#
#     em_out2 <- mixtools::multmixEM(m, theta = init_theta, lambda = init_lambda, k = 2)
#     list(k1 = em_out1, k2 = em_out2)
# }

fit_bm_cpp = function(x, n, ave_p, p1, p2, theta1, max_iter = 100, tol = 1e-6) {
    em_out = em_bm(x, n, p1, p2, theta1, max_iter, tol)

    loglike_k1 = binomial_log_likelihood(x, n, ave_p)
    pval_k1 = pbinom(x, n, ave_p, lower.tail = T)

    list(
        k1 = list(
            loglik = loglike_k1, 
            pval = pval_k1,
            p = ave_p 
        ),
        k2 = em_out
    )
}

#######################################################################
#                   End of internal function region                   #
#######################################################################


#' Fit binomial mixture distribution for one locus
#'
#' @param d_select_maj_base data.frame of one locus.
#' @param theta1 initial theta for the second component.
#' @param max_iter maximum iteration.
#' @param tol tolerance of log-likelihood to stop iteration.
#' @return list of p-value and model parameters.
#' @export
#' @examples
#'
#' ###
process_locus_bm = function(
    d_select_maj_base,
    theta1 = 0.9,
    max_iter = 100,
    tol = 1e-6
    ) {

    #################################################
    ### Transform data
    ## coverage
    n <- d_select_maj_base$coverage
    ## majority base depth
    x <- d_select_maj_base[, alt_depth]
    ## majority base forward depth
    y_fwd <- d_select_maj_base[, fwd_depth]
    ## majority base reverse depth
    y_rev <- d_select_maj_base[, rev_depth]
    ## average af
    mean_af <- sum(x) / sum(n)

    #################################################
    ### Model fitting
    # em_out_l = fit_bm_r(x, n, init_theta = NULL, init_lambda = NULL)
    # major_i = which.max(em_out_l$k2$lambda)
    # minor_i = which.min(em_out_l$k2$lambda)

    # list(
        # pval = em_out_l$k2$posterior[, major_i],
        # parameters = data.table(
            # loglik_k1 = em_out_l$k1$loglik,
            # loglik_k2 = em_out_l$k2$loglik,
            # k2_pi1 = em_out_l$k2$theta[major_i, 1],
            # k2_pi2 = em_out_l$k2$theta[minor_i, 1],
            # k2_theta1 = em_out_l$k2$lambda[major_i]
        # )
    # )

    ## use the binomial estimate as the initial value
    ave_p = sum(x)/sum(n)
    fit_l = fit_bm_cpp(x, n, ave_p = ave_p, p1 = ave_p, p2 = ave_p/1.1, theta1 = theta1, max_iter = max_iter, tol = tol)

    list(
        bm_pval = fit_l$k2$pval,
        bi_pval = fit_l$k1$pval,
        parameters = data.table(
            loglik_k1 = fit_l$k1$loglik,
            loglik_k2 = fit_l$k2$loglik,
            k1_pi = fit_l$k1$p,
            k2_pi1 = fit_l$k2$p1,
            k2_pi2 = fit_l$k2$p2,
            k2_theta1 = fit_l$k2$theta1
        )
    )
}

