#include <cmath>
#include <vector>
#include <Rcpp.h>
#include <Rmath.h>

using namespace std;

// Compute the log-likelihood of a binomial mixture model with given parameters
// [[Rcpp::export]]
double binomial_mixture_log_likelihood(
    const Rcpp::NumericVector x, 
    const Rcpp::NumericVector n, 
    double p1, 
    double p2, 
    double theta1
) {
    int N = x.size();
    double ll = 0.0;
    for (int i = 0; i < N; i++) {
        ll += log(theta1 * R::dbinom(x[i], n[i], p1, 0) + (1 - theta1) * R::dbinom(x[i], n[i], p2, 0));
    }
    return ll;
}

// Compute the log-likelihood of binomial model
// [[Rcpp::export]]
double binomial_log_likelihood(
    const Rcpp::NumericVector x,
    const Rcpp::NumericVector n,
    double p
) {
    int N = x.size();
    double ll = 0.0;
    for (int i = 0; i < N; i++) {
        ll += log(R::dbinom(x[i], n[i], p, 0));
    }
    return ll;
}


Rcpp::NumericVector estimate_posterior_p(
    const Rcpp::NumericVector x,
    const Rcpp::NumericVector n,
    double p
) {
    Rcpp::NumericVector pval = Rcpp::NumericVector(x.size());
    for (int i = 0; i < x.size(); i++) {
        pval[i] = R::pbinom(x[i], n[i], p, 1, 0);
    }
    return pval;
}

// Estimate the parameters of a binomial mixture model using the EM algorithm
// [[Rcpp::export]]
Rcpp::List em_bm(
    const Rcpp::NumericVector x,
    const Rcpp::NumericVector n,
    double p1 = 0.5, 
    double p2 = 0.6,
    double theta1 = 0.7,
    int max_iter = 1000, 
    double tol = 0.001
) {
    Rcpp::checkUserInterrupt();

    double N = n.size();
    // Run the EM algorithm

    double p1_old = 0.0, p2_old = 0.0;
    // parts for updating pi
    // parts for updating theta
    double theta_part1_deno = 0.0, theta_part1_num = 0.0; 
    double theta_part2_deno = 0.0, theta_part2_num = 0.0; 

    double ll = binomial_mixture_log_likelihood(x, n, p1, p2, theta1);
    double ll_old = 0.0;
    double pi_part1 = 0.0; 

    // Equation 6, term 2, and term 3
    // double term2_1_2, term3_1_2;

    // value of Equation 4 when k is 1 or 2
    double q1, q2;
    int iter = 0;

    // for (int iter = 0; iter < max_iter; iter++) {
    do {
        pi_part1 = 0.0; 
        p1_old = p1, p2_old = p2;
        ll_old = ll;


        for (int i = 0; i < N; i++) {
            // E-step: compute the responsibilities

            // term2_1_2 = pow(p2_old / p1_old, x[1]);
            // term3_1_2 = pow((1 - p2_old) / (1 - p1_old), n[i] - x[i]);
            // q1 = 1 / ((1-w) / w  * term2_1_2 * term3_1_2 + 1);
            // q2 = 1 / (w / (1 - w) * (1/term2_1_2) * (1/term3_1_2) + 1);

            q1 = 1 / ( (1-theta1) / theta1 * pow(p2_old / p1_old, x[i]) * pow((1 - p2_old) / (1 - p1_old), n[i] - x[i]) + 1);
            q2 = 1-q1;

            pi_part1 += q1;

            theta_part1_num += q1 * x[i];
            theta_part2_num += q2 * x[i];
            theta_part1_deno += q1 * n[i];
            theta_part2_deno += q2 * n[i];
        }

        // M-step: update the parameters
        theta1 = pi_part1 / N;
        p1 = theta_part1_num / theta_part1_deno;
        p2 = theta_part2_num / theta_part2_deno;

        // Rcpp::Rcout << "iter: " << iter << ", ll: " << ll << ", p1: " << pi_part1 << ", p2: " << pi_part2 << ", w: " << w << std::endl;
        // Rcpp::Rcout << "iter: " << iter << ", ll: " << ll << ", p1: " << p1 << ", p2: " << p2 << ", w: " << w << std::endl;

        // Check for convergence
        iter++;
        if (iter % 10 == 0 || iter < 10) {
            ll = binomial_mixture_log_likelihood(x, n, p1, p2, theta1);
            if (ll - ll_old < tol) {
                break;
            }
        }
    } while (iter < max_iter);

    Rcpp::NumericVector pval; 
    if (p1 < p2) {
        pval = estimate_posterior_p(x, n, p2);
    } else {
        pval = estimate_posterior_p(x, n, p1);
    }

    return Rcpp::List::create(
        Rcpp::Named("p1") = p1,
        Rcpp::Named("p2") = p2,
        Rcpp::Named("pval") = pval,
        Rcpp::Named("theta1") = theta1,
        Rcpp::Named("loglik") = ll,
        Rcpp::Named("iter") = iter
    );
}

