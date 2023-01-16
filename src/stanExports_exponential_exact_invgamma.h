// Generated by rstantools.  Do not edit by hand.

/*
    biPOD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    biPOD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with biPOD.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_exponential_exact_invgamma_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_exponential_exact_invgamma");
    reader.add_event(104, 102, "end", "model_exponential_exact_invgamma");
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
            stan::math::assign(prob, (prob + (((((stan::math::exp(binomial_coefficient_log(n, j)) * pow((1 - alpha_t), j)) * pow(alpha_t, (n - j))) * stan::math::exp(binomial_coefficient_log((m - 1), (m - j)))) * pow((1 - beta_t), j)) * pow(beta_t, (m - j)))));
        }
        current_statement_begin__ = 25;
        stan::math::assign(lprob, stan::math::log(prob));
        current_statement_begin__ = 26;
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
        current_statement_begin__ = 30;
        local_scalar_t__ prob(DUMMY_VAR__);
        (void) prob;  // dummy to suppress unused var warning
        stan::math::initialize(prob, DUMMY_VAR__);
        stan::math::fill(prob, DUMMY_VAR__);
        current_statement_begin__ = 31;
        local_scalar_t__ lprob(DUMMY_VAR__);
        (void) lprob;  // dummy to suppress unused var warning
        stan::math::initialize(lprob, DUMMY_VAR__);
        stan::math::fill(lprob, DUMMY_VAR__);
        current_statement_begin__ = 32;
        local_scalar_t__ a_t(DUMMY_VAR__);
        (void) a_t;  // dummy to suppress unused var warning
        stan::math::initialize(a_t, DUMMY_VAR__);
        stan::math::fill(a_t, DUMMY_VAR__);
        current_statement_begin__ = 33;
        local_scalar_t__ b_t(DUMMY_VAR__);
        (void) b_t;  // dummy to suppress unused var warning
        stan::math::initialize(b_t, DUMMY_VAR__);
        stan::math::fill(b_t, DUMMY_VAR__);
        current_statement_begin__ = 34;
        int upper_v(0);
        (void) upper_v;  // dummy to suppress unused var warning
        stan::math::fill(upper_v, std::numeric_limits<int>::min());
        current_statement_begin__ = 35;
        int lower_v(0);
        (void) lower_v;  // dummy to suppress unused var warning
        stan::math::fill(lower_v, std::numeric_limits<int>::min());
        current_statement_begin__ = 36;
        local_scalar_t__ x(DUMMY_VAR__);
        (void) x;  // dummy to suppress unused var warning
        stan::math::initialize(x, DUMMY_VAR__);
        stan::math::fill(x, DUMMY_VAR__);
        current_statement_begin__ = 37;
        local_scalar_t__ cdf(DUMMY_VAR__);
        (void) cdf;  // dummy to suppress unused var warning
        stan::math::initialize(cdf, DUMMY_VAR__);
        stan::math::fill(cdf, DUMMY_VAR__);
        current_statement_begin__ = 38;
        int m(0);
        (void) m;  // dummy to suppress unused var warning
        stan::math::fill(m, std::numeric_limits<int>::min());
        current_statement_begin__ = 40;
        stan::math::assign(x, uniform_rng(0, 1, base_rng__));
        current_statement_begin__ = 41;
        stan::math::assign(cdf, 0);
        current_statement_begin__ = 42;
        stan::math::assign(m, -(1));
        current_statement_begin__ = 44;
        while (as_bool(logical_lt(cdf, x))) {
            current_statement_begin__ = 45;
            stan::math::assign(m, (m + 1));
            current_statement_begin__ = 47;
            stan::math::assign(upper_v, (n - 1));
            current_statement_begin__ = 48;
            if (as_bool(logical_gte((n - m), 0))) {
                current_statement_begin__ = 49;
                stan::math::assign(lower_v, (n - m));
            } else {
                current_statement_begin__ = 51;
                stan::math::assign(lower_v, 0);
            }
            current_statement_begin__ = 54;
            stan::math::assign(a_t, ((mu * (stan::math::exp((delta_t * (lambda - mu))) - 1)) / ((lambda * stan::math::exp(((lambda - mu) * delta_t))) - mu)));
            current_statement_begin__ = 55;
            stan::math::assign(b_t, ((lambda * a_t) / mu));
            current_statement_begin__ = 57;
            stan::math::assign(prob, 0);
            current_statement_begin__ = 59;
            for (int j = lower_v; j <= upper_v; ++j) {
                current_statement_begin__ = 60;
                stan::math::assign(prob, (prob + ((((stan::math::exp(binomial_coefficient_log(n, j)) * stan::math::exp(binomial_coefficient_log((m - 1), ((n - j) - 1)))) * pow(a_t, j)) * pow(((1 - a_t) * (1 - b_t)), (n - j))) * pow(b_t, ((m - n) + j)))));
            }
            current_statement_begin__ = 63;
            stan::math::assign(cdf, (cdf + prob));
        }
        current_statement_begin__ = 66;
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
#include <stan_meta_header.hpp>
class model_exponential_exact_invgamma
  : public stan::model::model_base_crtp<model_exponential_exact_invgamma> {
private:
        int S;
        int n0;
        double t0;
        std::vector<int> N;
        std::vector<double> T;
        double a;
        double b;
public:
    model_exponential_exact_invgamma(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_exponential_exact_invgamma(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_exponential_exact_invgamma_namespace::model_exponential_exact_invgamma";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 71;
            context__.validate_dims("data initialization", "S", "int", context__.to_vec());
            S = int(0);
            vals_i__ = context__.vals_i("S");
            pos__ = 0;
            S = vals_i__[pos__++];
            check_greater_or_equal(function__, "S", S, 2);
            current_statement_begin__ = 72;
            context__.validate_dims("data initialization", "n0", "int", context__.to_vec());
            n0 = int(0);
            vals_i__ = context__.vals_i("n0");
            pos__ = 0;
            n0 = vals_i__[pos__++];
            check_greater_or_equal(function__, "n0", n0, 0);
            current_statement_begin__ = 73;
            context__.validate_dims("data initialization", "t0", "double", context__.to_vec());
            t0 = double(0);
            vals_r__ = context__.vals_r("t0");
            pos__ = 0;
            t0 = vals_r__[pos__++];
            check_greater_or_equal(function__, "t0", t0, 0);
            current_statement_begin__ = 75;
            validate_non_negative_index("N", "S", S);
            context__.validate_dims("data initialization", "N", "int", context__.to_vec(S));
            N = std::vector<int>(S, int(0));
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            size_t N_k_0_max__ = S;
            for (size_t k_0__ = 0; k_0__ < N_k_0_max__; ++k_0__) {
                N[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 76;
            validate_non_negative_index("T", "S", S);
            context__.validate_dims("data initialization", "T", "double", context__.to_vec(S));
            T = std::vector<double>(S, double(0));
            vals_r__ = context__.vals_r("T");
            pos__ = 0;
            size_t T_k_0_max__ = S;
            for (size_t k_0__ = 0; k_0__ < T_k_0_max__; ++k_0__) {
                T[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 78;
            context__.validate_dims("data initialization", "a", "double", context__.to_vec());
            a = double(0);
            vals_r__ = context__.vals_r("a");
            pos__ = 0;
            a = vals_r__[pos__++];
            check_greater_or_equal(function__, "a", a, 0);
            current_statement_begin__ = 79;
            context__.validate_dims("data initialization", "b", "double", context__.to_vec());
            b = double(0);
            vals_r__ = context__.vals_r("b");
            pos__ = 0;
            b = vals_r__[pos__++];
            check_greater_or_equal(function__, "b", b, 0);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 83;
            num_params_r__ += 1;
            current_statement_begin__ = 84;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_exponential_exact_invgamma() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 83;
        if (!(context__.contains_r("lambda")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable lambda missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("lambda");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "lambda", "double", context__.to_vec());
        double lambda(0);
        lambda = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, lambda);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable lambda: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 84;
        if (!(context__.contains_r("mu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable mu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("mu");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "mu", "double", context__.to_vec());
        double mu(0);
        mu = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, mu);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable mu: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 83;
            local_scalar_t__ lambda;
            (void) lambda;  // dummy to suppress unused var warning
            if (jacobian__)
                lambda = in__.scalar_lb_constrain(0, lp__);
            else
                lambda = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 84;
            local_scalar_t__ mu;
            (void) mu;  // dummy to suppress unused var warning
            if (jacobian__)
                mu = in__.scalar_lb_constrain(0, lp__);
            else
                mu = in__.scalar_lb_constrain(0);
            // transformed parameters
            current_statement_begin__ = 88;
            local_scalar_t__ ro;
            (void) ro;  // dummy to suppress unused var warning
            stan::math::initialize(ro, DUMMY_VAR__);
            stan::math::fill(ro, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 89;
            stan::math::assign(ro, (lambda - mu));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 88;
            if (stan::math::is_uninitialized(ro)) {
                std::stringstream msg__;
                msg__ << "Undefined transformed parameter: ro";
                stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable ro: ") + msg__.str()), current_statement_begin__, prog_reader__());
            }
            // model body
            current_statement_begin__ = 93;
            lp_accum__.add(inv_gamma_log<propto__>(lambda, a, b));
            current_statement_begin__ = 94;
            lp_accum__.add(inv_gamma_log<propto__>(mu, a, b));
            current_statement_begin__ = 96;
            for (int i = 1; i <= S; ++i) {
                current_statement_begin__ = 97;
                lp_accum__.add(birthDeathLike_log<propto__>(get_base1(N, i, "N", 1), n0, (get_base1(T, i, "T", 1) - t0), lambda, mu, pstream__));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("lambda");
        names__.push_back("mu");
        names__.push_back("ro");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_exponential_exact_invgamma_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double lambda = in__.scalar_lb_constrain(0);
        vars__.push_back(lambda);
        double mu = in__.scalar_lb_constrain(0);
        vars__.push_back(mu);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 88;
            double ro;
            (void) ro;  // dummy to suppress unused var warning
            stan::math::initialize(ro, DUMMY_VAR__);
            stan::math::fill(ro, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 89;
            stan::math::assign(ro, (lambda - mu));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                vars__.push_back(ro);
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_exponential_exact_invgamma";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "lambda";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "ro";
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "lambda";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "ro";
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_exponential_exact_invgamma_namespace::model_exponential_exact_invgamma stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
