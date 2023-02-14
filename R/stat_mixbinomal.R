library(mixtools)

fit_bm <- function(y, N, init_theta, init_lambda, max_try = 1) {

    m = data.matrix(data.frame(y, N-y))
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



process_locus_bm = function(d_select_maj_base, max_try = 1) {

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
    ### Model fitting
    ## TODO: Set init value
    em_out_l = fit_bm(y, N, init_theta = matrix(rep(0.5, 4), 2, 2), init_lambda = c(0.2, 0.8), max_try = max_try)

    pval_m = em_out_l$k2$posterior
    #     d_select_maj_base$p1 = pval_m[, 1]
    #     d_select_maj_base$p2 = pval_m[, 2]
    ## choose the p with bigger lambda/pi, which means the p for H0
    pval = pval_m[, which.max(em_out_l$k2$lambda)]

    list(
        #         res = d_select_maj_base,
        pval = pval,
        parameters = data.table(
            loglik_k1 = em_out_l$k1$loglik,
            loglik_k2 = em_out_l$k2$loglik,
            k2_pi1 = em_out_l$k2$lambda[1],
            k2_pi2 = em_out_l$k2$lambda[2],
            k2_theta1 = em_out_l$k2$theta[1, 1],
            k2_theta2 = em_out_l$k2$theta[2, 1]
        )
    )
}

