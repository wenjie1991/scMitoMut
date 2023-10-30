#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;

// calculate parital derivative of log-likelihood function with respect to a
double ll_a(
    const double a, const double b,
    const long * x, const long * n,
    const long N
) {
    double term1 = N * (R::digamma(a + b) - R::digamma(a));

    double term2 = 0;
    for (int i = 0; i < N; i++) {
        term2 += R::digamma(a + b + n[i]);
    }

    double term3 = 0;
    for (int i = 0; i < N; i++) {
        // term3 += R::digamma(a + x[i]);
        term3 += R::digamma(a + x[i]);
    }

    return term1 - term2 + term3;
}

// calculate parital derivative of log-likelihood function with respect to b
double ll_b(
    const double a, const double b,
    const long * x, const long * n,
    const long N
) {
    double term1 = N * (R::digamma(a + b) - R::digamma(b));

    double term2 = 0;
    for (int i = 0; i < N; i++) {
        term2 += R::digamma(a + b + n[i]);
    }

    double term3 = 0;
    for (int i = 0; i < N; i++) {
        term3 += R::digamma(b + n[i] - x[i]);
    }

    return term1 - term2 + term3;
}

// calculate paritial partial derivative of log-likelihood function with respect to a
double ll_aa(
    const double a, const double b,
    const long * x, const long * n,
    const long N
) {
    double term1 = N * (R::trigamma(a + b) - R::trigamma(a));

    double term2 = 0;
    for (int i = 0; i < N; i++) {
        term2 += R::trigamma(a + b + n[i]);
    }

    double term3 = 0;
    for (int i = 0; i < N; i++) {
        term3 += R::trigamma(a + x[i]);
    }

    return term1 - term2 + term3;
}

// calculate paritial partial derivative of log-likelihood function with respect to b
double ll_bb(
    const double a, const double b,
    const long * x, const long * n,
    const long N
) {
    double term1 = N * (R::trigamma(a + b) - R::trigamma(b));

    double term2 = 0;
    for (int i = 0; i < N; i++) {
        term2 += R::trigamma(a + b + n[i]);
    }

    double term3 = 0;
    for (int i = 0; i < N; i++) {
        term3 += R::trigamma(b + n[i] - x[i]);
    }

    return term1 - term2 + term3;
}

// calculate paritial partial derivative of log-likelihood function with respect to a and b
double ll_ab(
    const double a, const double b,
    const long * x, const long * n,
    const long N
) {
    double term1 = N * R::trigamma(a + b);

    double term2 = 0;
    for (int i = 0; i < N; i++) {
        term2 += R::trigamma(a + b + n[i]);
    }

    return term1 - term2;
}

// 1st moment estimation of beta-binomial distribution
double m1( const long * x, const long N) {
    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += x[i];
    }
    return sum / N;
}

// 2nd moment estimation of beta-binomial distribution
double m2(const long * x, const long N) {
    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += pow(x[i], 2);
    }
    return sum / N;
}

// moment estimation of beta-binomial distribution and return a
double m_a(const double m1, const double m2, const double n_mean) {
    return abs(m1 * (n_mean * m1 - m2) / (n_mean * (m2 - pow(m1, 2) - m1) + pow(m1, 2)));
}

// moment estimation of beta-binomial distribution and return b
double m_b(const double m1, const double m2, const double n_mean) {
    return abs((n_mean - m1) * (n_mean * m1 - m2) / (n_mean * (m2 - pow(m1,2) - m1) + pow(m1,2)));
}

// log n choose k
// [[Rcpp::export]]
double log_n_choose_k(long n, long k) {
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 0;

    double result = log(n);
    for(long i = 2; i <= k; ++i ) {
        result += log(n-i+1);
        result -= log(i);
    }
    return result;
}

long double log_combination(int n, int k) {
    return R::lgammafn(n + 1) - R::lgammafn(k + 1) - R::lgammafn(n - k + 1);
}

long double log_beta(double a, double b) {
    return R::lgammafn(a) + R::lgammafn(b) - R::lgammafn(a + b);
}

// likelihood beta-binomial distribution
// [[Rcpp::export]]
long double log_beta_binomial_pmf(
    const double a, const double b,
    const double x, const double n
) {
    // return log_n_choose_k(n, x) + log(R::beta(x + a, n - x + b)) -  log(R::beta(a, b));
    return log_combination(n, x) + log_beta(x + a, n - x + b) -  log_beta(a, b);
}


long double pbetabinom_c(
	double x, const double n,
    const double a, const double b
) {
    long double p = 0;
    if (2 * x > n) {
        for (int j = x + 1; j <= n; j++) {
            p += exp(log_beta_binomial_pmf(a, b, j, n));
        }
        p = 1-p;
    } else {
        for (int j = 0; j <= x; j++) {
            p += exp(log_beta_binomial_pmf(a, b, j, n));
        }
    }
    return p;
}
	

// Calculate p value for beta-binomial distribution
// [[Rcpp::export]]
Rcpp::NumericVector pbetabinom(
    const Rcpp::NumericVector x, const Rcpp::NumericVector n,
    const double a, const double b
) {
    int N = x.size();
    Rcpp::NumericVector p(N);
    for (int i = 0; i < N; i++) {
        p[i] = pbetabinom_c(x[i], n[i], a, b);
    }
    return p;
}

// log-likelihood function
double LL_beta(
    const double a, const double b,
    const long * x, const long * n,
    const long N
) {
    double ll = 0;
    for (int i = 0; i < N; i++) {
        // Debug beta-binomial distribution
        // double tmp = log_beta_binomial_pmf(a, b, x[i], n[i]);
        // if (tmp > 0) {
            // Rcpp::Rcout << "a: " << a << " b: " << b << " x: " << x[i] << " n: " << n[i] << " N: " << N << std::endl;
        // }
        ll += log_beta_binomial_pmf(a, b, x[i], n[i]);
    }
    return ll;
}


// mle beta-binomial distribution
// [[Rcpp::export]]
Rcpp::List mle_bb(
    Rcpp::IntegerVector x, Rcpp::IntegerVector n,
    int max_iter = 100, double tol = 1e-3
) {
    Rcpp::checkUserInterrupt();

    long * x_ptr = new long[x.size()];
    long * n_ptr = new long[n.size()];
    long N = x.size();

    for (int i = 0; i < N; i++) {
        x_ptr[i] = x[i];
        n_ptr[i] = n[i];
    }

    double mean = m1(x_ptr, N);
    double var = m2(x_ptr, N);
    double n_mean = 0;
    for (int i = 0; i < N; i++) {
        n_mean += n_ptr[i];
    }
    n_mean /= N;

    // moment estimation
    double a = m_a(mean, var, n_mean);
    double b = m_b(mean, var, n_mean);

    // newton-raphson method
    double a_new = a, b_new = b;
    double a_old, b_old;
    double ll_old, ll_new;
    ll_new = LL_beta(a_new, b_new, x_ptr, n_ptr, N);
    int iter = 0;
    do {
        a_old = a_new;
        b_old = b_new;

        ll_old = ll_new;

        // parameter vector
        arma::mat parm(2, 1);
        parm(0, 0) = a_old;
        parm(1, 0) = b_old;

        // partial derivative vector
        arma::mat partial_v(2, 1);
        partial_v(0, 0) = ll_a(a_old, b_old, x_ptr, n_ptr, N);
        partial_v(1, 0) = ll_b(a_old, b_old, x_ptr, n_ptr, N);

        // Jacobian matrix
        arma::mat J_pre(2, 2);
        J_pre(0, 0) = ll_aa(a_old, b_old, x_ptr, n_ptr, N);
        J_pre(0, 1) = ll_ab(a_old, b_old, x_ptr, n_ptr, N);
        J_pre(1, 0) = ll_ab(a_old, b_old, x_ptr, n_ptr, N);
        J_pre(1, 1) = ll_bb(a_old, b_old, x_ptr, n_ptr, N);
        arma::mat J = arma::inv(J_pre);

        // update a and b
        arma::mat f2 = parm - J * partial_v;

        a_new = f2(0, 0);
        b_new = f2(1, 0);

        if (a_new < 0) a_new = a;
        if (b_new < 0) b_new = b;

        // std::cout << ll_new << std::endl;

        ll_new = LL_beta(a_new, b_new, x_ptr, n_ptr, N);

        iter++;
    // } while (ll_new - ll_old > tol && ll_new < -tol && iter < max_iter);
    } while (ll_new - ll_old > tol && iter < max_iter);
    // } while (abs(_new - a_old) > tol && abs(b_new - b_old) > tol && iter < max_iter);

    return Rcpp::List::create(
        Rcpp::Named("a") = a_new,
        Rcpp::Named("b") = b_new,
        Rcpp::Named("mean") = mean,
        Rcpp::Named("var") = var,
        Rcpp::Named("n_mean") = n_mean,
        Rcpp::Named("iter") = iter
    );
}
