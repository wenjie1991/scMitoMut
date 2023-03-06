library(mixtools)

fit_bm_r <- function(x, n, init_theta, init_lambda, max_try = 1) {

    m = data.matrix(data.frame(x, n-x))
    em_out1 <- mixtools::multmixEM(m, k = 1)

    try_num <- 0
    while (try_num < max_try) {
        try_num <- try_num + 1
        tryCatch(
            {
                em_out2 <- mixtools::multmixEM(m, theta = init_theta, lambda = init_lambda, k = 2)
                break
            },
            error = function(e) {
                print(e)
                print(paste0("Try ", try_num, " times."))
            }
        )
        init_vaf <- init_vaf - 0.001
    }
    list(k1 = em_out1, k2 = em_out2)
}

fit_bm_cpp = function(x, n, ave_p, p1, p2, theta1, max_iter = 100, tol = 1e-6) {
    em_out = em_bm(x, n, p1, p2, theta1, max_iter, tol)

    loglike_k1 = binomial_log_likelihood(x, n, ave_p)

    list(k1_ll = loglike_k1, k2 = em_out)
}



process_locus_bm = function(d_select_maj_base, max_try = 1) {

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
    ## averrage vaf
    mean_vaf <- sum(x) / sum(n)

    #################################################
    ### Model fitting
    # em_out_l = fit_bm_r(x, n, init_theta = matrix(rep(0.5, 4), 2, 2), init_lambda = c(0.2, 0.8), max_try = max_try)
    # pval_m = em_out_l$k2$posterior
    # pval = pval_m[, which.max(em_out_l$k2$lambda)]
    # list(
    #     #         res = d_select_maj_base,
    #     pval = pval,
    #     parameters = data.table(
    #         loglik_k1 = em_out_l$k1$loglik,
    #         loglik_k2 = em_out_l$k2$loglik,
    #         k2_pi1 = em_out_l$k2$lambda[1],
    #         k2_pi2 = em_out_l$k2$lambda[2],
    #         k2_theta1 = em_out_l$k2$theta[1, 1],
    #         k2_theta2 = em_out_l$k2$theta[2, 1]
    #     )
    # )

    ## use the binomial estimate as the initial value
    ave_p = sum(x)/sum(n)
    fit_l = fit_bm_cpp(x, n, ave_p = ave_p, p1 = ave_p, p2 = ave_p/2, theta1 = 0.95, max_iter = 100, tol = 1e-6)

    list(
        pval = fit_l$k2$pval,
        parameters = data.table(
            loglik_k1 = fit_l$k1_ll,
            loglik_k2 = fit_l$k2$loglik,
            k2_pi1 = fit_l$k2$p1,
            k2_pi2 = fit_l$k2$p2,
            k2_theta1 = fit_l$k2$theta1
        )
    )
}

