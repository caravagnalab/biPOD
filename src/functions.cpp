// Generated by rstantools.  Do not edit by hand.

// [[Rcpp::depends(StanHeaders)]]
// [[Rcpp::depends(rstan)]]
// [[Rcpp::plugins(rstan)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <boost/integer/integer_log2.hpp>
#include <RcppEigen.h>
// Code generated by Stan version 2.21.0
#include <stan/model/standalone_functions_header.hpp>
namespace user_be3f9d9c2623925b39a604ca5a538e65 { 
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using namespace stan::math;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_d;
typedef Eigen::Matrix<double, 1, Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_d;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "unknown file name");
    reader.add_event(72, 70, "end", "unknown file name");
    return reader;
}
template <bool propto, typename T2__, typename T3__, typename T4__>
typename boost::math::tools::promote_args<T2__, T3__, T4__>::type
birthDeathLike_log(const int& m,
                       const int& n,
                       const T2__& delta_t,
                       const T3__& lambda,
                       const T4__& mu, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T2__, T3__, T4__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 4;
        local_scalar_t__ prob(DUMMY_VAR__);
        (void) prob;  // dummy to suppress unused var warning
        stan::math::initialize(prob, DUMMY_VAR__);
        stan::math::fill(prob, DUMMY_VAR__);
        current_statement_begin__ = 5;
        local_scalar_t__ lprob(DUMMY_VAR__);
        (void) lprob;  // dummy to suppress unused var warning
        stan::math::initialize(lprob, DUMMY_VAR__);
        stan::math::fill(lprob, DUMMY_VAR__);
        current_statement_begin__ = 6;
        local_scalar_t__ alpha_t(DUMMY_VAR__);
        (void) alpha_t;  // dummy to suppress unused var warning
        stan::math::initialize(alpha_t, DUMMY_VAR__);
        stan::math::fill(alpha_t, DUMMY_VAR__);
        current_statement_begin__ = 7;
        local_scalar_t__ beta_t(DUMMY_VAR__);
        (void) beta_t;  // dummy to suppress unused var warning
        stan::math::initialize(beta_t, DUMMY_VAR__);
        stan::math::fill(beta_t, DUMMY_VAR__);
        current_statement_begin__ = 8;
        int upper_sum_value(0);
        (void) upper_sum_value;  // dummy to suppress unused var warning
        stan::math::fill(upper_sum_value, std::numeric_limits<int>::min());
        current_statement_begin__ = 10;
        if (as_bool(logical_lte(n, m))) {
            current_statement_begin__ = 11;
            stan::math::assign(upper_sum_value, n);
        } else {
            current_statement_begin__ = 13;
            stan::math::assign(upper_sum_value, m);
        }
        current_statement_begin__ = 16;
        stan::math::assign(alpha_t, ((mu * (stan::math::exp((delta_t * (lambda - mu))) - 1)) / ((lambda * stan::math::exp(((lambda - mu) * delta_t))) - mu)));
        current_statement_begin__ = 17;
        stan::math::assign(beta_t, ((lambda * alpha_t) / mu));
        current_statement_begin__ = 19;
        stan::math::assign(prob, 0);
        current_statement_begin__ = 21;
        for (int j = 1; j <= upper_sum_value; ++j) {
            current_statement_begin__ = 22;
            if (as_bool((primitive_value(logical_neq(stan::math::exp(binomial_coefficient_log(n, j)), stan::math::positive_infinity())) && primitive_value(logical_neq(stan::math::exp(binomial_coefficient_log((m - 1), (m - j))), stan::math::positive_infinity()))))) {
                current_statement_begin__ = 23;
                stan::math::assign(prob, (prob + (((((stan::math::exp(binomial_coefficient_log(n, j)) * pow((1 - alpha_t), j)) * pow(alpha_t, (n - j))) * stan::math::exp(binomial_coefficient_log((m - 1), (m - j)))) * pow((1 - beta_t), j)) * pow(beta_t, (m - j)))));
            }
        }
        current_statement_begin__ = 27;
        stan::math::assign(lprob, stan::math::log(prob));
        current_statement_begin__ = 28;
        return stan::math::promote_scalar<fun_return_scalar_t__>(lprob);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
template <typename T2__, typename T3__, typename T4__>
typename boost::math::tools::promote_args<T2__, T3__, T4__>::type
birthDeathLike_log(const int& m,
                       const int& n,
                       const T2__& delta_t,
                       const T3__& lambda,
                       const T4__& mu, std::ostream* pstream__) {
    return birthDeathLike_log<false>(m,n,delta_t,lambda,mu, pstream__);
}
struct birthDeathLike_log_functor__ {
    template <bool propto, typename T2__, typename T3__, typename T4__>
        typename boost::math::tools::promote_args<T2__, T3__, T4__>::type
    operator()(const int& m,
                       const int& n,
                       const T2__& delta_t,
                       const T3__& lambda,
                       const T4__& mu, std::ostream* pstream__) const {
        return birthDeathLike_log(m, n, delta_t, lambda, mu, pstream__);
    }
};
template <typename T1__, typename T2__, typename T3__, class RNG>
int
birthDeathLike_rng(const int& n,
                       const T1__& delta_t,
                       const T2__& lambda,
                       const T3__& mu, RNG& base_rng__, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T1__, T2__, T3__>::type local_scalar_t__;
    typedef int fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 32;
        local_scalar_t__ prob(DUMMY_VAR__);
        (void) prob;  // dummy to suppress unused var warning
        stan::math::initialize(prob, DUMMY_VAR__);
        stan::math::fill(prob, DUMMY_VAR__);
        current_statement_begin__ = 33;
        local_scalar_t__ lprob(DUMMY_VAR__);
        (void) lprob;  // dummy to suppress unused var warning
        stan::math::initialize(lprob, DUMMY_VAR__);
        stan::math::fill(lprob, DUMMY_VAR__);
        current_statement_begin__ = 34;
        local_scalar_t__ a_t(DUMMY_VAR__);
        (void) a_t;  // dummy to suppress unused var warning
        stan::math::initialize(a_t, DUMMY_VAR__);
        stan::math::fill(a_t, DUMMY_VAR__);
        current_statement_begin__ = 35;
        local_scalar_t__ b_t(DUMMY_VAR__);
        (void) b_t;  // dummy to suppress unused var warning
        stan::math::initialize(b_t, DUMMY_VAR__);
        stan::math::fill(b_t, DUMMY_VAR__);
        current_statement_begin__ = 36;
        int upper_v(0);
        (void) upper_v;  // dummy to suppress unused var warning
        stan::math::fill(upper_v, std::numeric_limits<int>::min());
        current_statement_begin__ = 37;
        int lower_v(0);
        (void) lower_v;  // dummy to suppress unused var warning
        stan::math::fill(lower_v, std::numeric_limits<int>::min());
        current_statement_begin__ = 38;
        local_scalar_t__ x(DUMMY_VAR__);
        (void) x;  // dummy to suppress unused var warning
        stan::math::initialize(x, DUMMY_VAR__);
        stan::math::fill(x, DUMMY_VAR__);
        current_statement_begin__ = 39;
        local_scalar_t__ cdf(DUMMY_VAR__);
        (void) cdf;  // dummy to suppress unused var warning
        stan::math::initialize(cdf, DUMMY_VAR__);
        stan::math::fill(cdf, DUMMY_VAR__);
        current_statement_begin__ = 40;
        int m(0);
        (void) m;  // dummy to suppress unused var warning
        stan::math::fill(m, std::numeric_limits<int>::min());
        current_statement_begin__ = 42;
        stan::math::assign(x, uniform_rng(0, 1, base_rng__));
        current_statement_begin__ = 43;
        stan::math::assign(cdf, 0);
        current_statement_begin__ = 44;
        stan::math::assign(m, -(1));
        current_statement_begin__ = 46;
        while (as_bool(logical_lt(cdf, x))) {
            current_statement_begin__ = 47;
            stan::math::assign(m, (m + 1));
            current_statement_begin__ = 49;
            stan::math::assign(upper_v, (n - 1));
            current_statement_begin__ = 50;
            if (as_bool(logical_gte((n - m), 0))) {
                current_statement_begin__ = 51;
                stan::math::assign(lower_v, (n - m));
            } else {
                current_statement_begin__ = 53;
                stan::math::assign(lower_v, 0);
            }
            current_statement_begin__ = 56;
            stan::math::assign(a_t, ((mu * (stan::math::exp((delta_t * (lambda - mu))) - 1)) / ((lambda * stan::math::exp(((lambda - mu) * delta_t))) - mu)));
            current_statement_begin__ = 57;
            stan::math::assign(b_t, ((lambda * a_t) / mu));
            current_statement_begin__ = 59;
            stan::math::assign(prob, 0);
            current_statement_begin__ = 61;
            for (int j = lower_v; j <= upper_v; ++j) {
                current_statement_begin__ = 62;
                stan::math::assign(prob, (prob + ((((stan::math::exp(binomial_coefficient_log(n, j)) * stan::math::exp(binomial_coefficient_log((m - 1), ((n - j) - 1)))) * pow(a_t, j)) * pow(((1 - a_t) * (1 - b_t)), (n - j))) * pow(b_t, ((m - n) + j)))));
            }
            current_statement_begin__ = 65;
            stan::math::assign(cdf, (cdf + prob));
        }
        current_statement_begin__ = 68;
        return stan::math::promote_scalar<fun_return_scalar_t__>(m);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct birthDeathLike_rng_functor__ {
    template <typename T1__, typename T2__, typename T3__, class RNG>
        int
    operator()(const int& n,
                       const T1__& delta_t,
                       const T2__& lambda,
                       const T3__& mu, RNG& base_rng__, std::ostream* pstream__) const {
        return birthDeathLike_rng(n, delta_t, lambda, mu, base_rng__, pstream__);
    }
};
 } 
// [[Rcpp::export]]
double
birthDeathLike_log(const int& m,
                       const int& n,
                       const double& delta_t,
                       const double& lambda,
                       const double& mu, std::ostream* pstream__ = 0){
  return 
user_be3f9d9c2623925b39a604ca5a538e65::birthDeathLike_log<false, double, double, double>(m, n, delta_t, lambda, mu, pstream__);
}
// [[Rcpp::export]]
int
birthDeathLike_rng(const int& n,
                       const double& delta_t,
                       const double& lambda,
                       const double& mu, boost::ecuyer1988& base_rng__, std::ostream* pstream__ = 0){
  return 
user_be3f9d9c2623925b39a604ca5a538e65::birthDeathLike_rng<double, double, double, boost::ecuyer1988>(n, delta_t, lambda, mu, base_rng__, pstream__);
}
