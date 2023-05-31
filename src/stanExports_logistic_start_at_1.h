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
namespace model_logistic_start_at_1_namespace {
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
    reader.add_event(0, 0, "start", "model_logistic_start_at_1");
    reader.add_event(82, 80, "end", "model_logistic_start_at_1");
    return reader;
}
template <typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
logistic_growth(const T0__& t,
                    const T1__& n0,
                    const T2__& rho,
                    const T3__& K, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 4;
        local_scalar_t__ num(DUMMY_VAR__);
        (void) num;  // dummy to suppress unused var warning
        stan::math::initialize(num, DUMMY_VAR__);
        stan::math::fill(num, DUMMY_VAR__);
        stan::math::assign(num,(K * n0));
        current_statement_begin__ = 5;
        local_scalar_t__ den(DUMMY_VAR__);
        (void) den;  // dummy to suppress unused var warning
        stan::math::initialize(den, DUMMY_VAR__);
        stan::math::fill(den, DUMMY_VAR__);
        stan::math::assign(den,(n0 + ((K - n0) * stan::math::exp((-(rho) * t)))));
        current_statement_begin__ = 6;
        return stan::math::promote_scalar<fun_return_scalar_t__>((num / den));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct logistic_growth_functor__ {
    template <typename T0__, typename T1__, typename T2__, typename T3__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
    operator()(const T0__& t,
                    const T1__& n0,
                    const T2__& rho,
                    const T3__& K, std::ostream* pstream__) const {
        return logistic_growth(t, n0, rho, K, pstream__);
    }
};
template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__, typename T5__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__>::type>::type
mean_t(const T0__& t,
           const T1__& t0,
           const T2__& n0,
           const std::vector<T3__>& t_array,
           const std::vector<T4__>& rho_array,
           const T5__& K, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__>::type>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 10;
        local_scalar_t__ current_n0(DUMMY_VAR__);
        (void) current_n0;  // dummy to suppress unused var warning
        stan::math::initialize(current_n0, DUMMY_VAR__);
        stan::math::fill(current_n0, DUMMY_VAR__);
        stan::math::assign(current_n0,n0);
        current_statement_begin__ = 11;
        local_scalar_t__ dt(DUMMY_VAR__);
        (void) dt;  // dummy to suppress unused var warning
        stan::math::initialize(dt, DUMMY_VAR__);
        stan::math::fill(dt, DUMMY_VAR__);
        current_statement_begin__ = 12;
        int n_t(0);
        (void) n_t;  // dummy to suppress unused var warning
        stan::math::fill(n_t, std::numeric_limits<int>::min());
        stan::math::assign(n_t,num_elements(t_array));
        current_statement_begin__ = 13;
        int n_rho(0);
        (void) n_rho;  // dummy to suppress unused var warning
        stan::math::fill(n_rho, std::numeric_limits<int>::min());
        stan::math::assign(n_rho,num_elements(rho_array));
        current_statement_begin__ = 15;
        if (as_bool(logical_eq(n_t, 0))) {
            current_statement_begin__ = 16;
            stan::math::assign(dt, (t - t0));
            current_statement_begin__ = 17;
            return stan::math::promote_scalar<fun_return_scalar_t__>(logistic_growth(dt, current_n0, get_base1(rho_array, 1, "rho_array", 1), K, pstream__));
        }
        current_statement_begin__ = 20;
        if (as_bool(logical_lte(t, get_base1(t_array, 1, "t_array", 1)))) {
            current_statement_begin__ = 21;
            stan::math::assign(dt, (t - t0));
            current_statement_begin__ = 22;
            return stan::math::promote_scalar<fun_return_scalar_t__>(logistic_growth(dt, current_n0, get_base1(rho_array, 1, "rho_array", 1), K, pstream__));
        }
        current_statement_begin__ = 25;
        stan::math::assign(dt, (get_base1(t_array, 1, "t_array", 1) - t0));
        current_statement_begin__ = 26;
        stan::math::assign(current_n0, logistic_growth(dt, current_n0, get_base1(rho_array, 1, "rho_array", 1), K, pstream__));
        current_statement_begin__ = 27;
        for (int i = 2; i <= n_t; ++i) {
            current_statement_begin__ = 28;
            if (as_bool(logical_lte(t, get_base1(t_array, i, "t_array", 1)))) {
                current_statement_begin__ = 29;
                stan::math::assign(dt, (t - get_base1(t_array, (i - 1), "t_array", 1)));
                current_statement_begin__ = 30;
                return stan::math::promote_scalar<fun_return_scalar_t__>(logistic_growth(dt, current_n0, get_base1(rho_array, i, "rho_array", 1), K, pstream__));
            } else {
                current_statement_begin__ = 32;
                stan::math::assign(dt, (get_base1(t_array, i, "t_array", 1) - get_base1(t_array, (i - 1), "t_array", 1)));
                current_statement_begin__ = 33;
                stan::math::assign(current_n0, logistic_growth(dt, current_n0, get_base1(rho_array, i, "rho_array", 1), K, pstream__));
            }
        }
        current_statement_begin__ = 37;
        stan::math::assign(dt, (t - get_base1(t_array, n_t, "t_array", 1)));
        current_statement_begin__ = 38;
        return stan::math::promote_scalar<fun_return_scalar_t__>(logistic_growth(dt, current_n0, get_base1(rho_array, n_rho, "rho_array", 1), K, pstream__));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct mean_t_functor__ {
    template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__, typename T5__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__>::type>::type
    operator()(const T0__& t,
           const T1__& t0,
           const T2__& n0,
           const std::vector<T3__>& t_array,
           const std::vector<T4__>& rho_array,
           const T5__& K, std::ostream* pstream__) const {
        return mean_t(t, t0, n0, t_array, rho_array, K, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_logistic_start_at_1
  : public stan::model::model_base_crtp<model_logistic_start_at_1> {
private:
        int S;
        int G;
        std::vector<int> N;
        std::vector<double> T;
        std::vector<double> t_array;
        double prior_K;
public:
    model_logistic_start_at_1(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_logistic_start_at_1(stan::io::var_context& context__,
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
        static const char* function__ = "model_logistic_start_at_1_namespace::model_logistic_start_at_1";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 43;
            context__.validate_dims("data initialization", "S", "int", context__.to_vec());
            S = int(0);
            vals_i__ = context__.vals_i("S");
            pos__ = 0;
            S = vals_i__[pos__++];
            check_greater_or_equal(function__, "S", S, 1);
            current_statement_begin__ = 44;
            context__.validate_dims("data initialization", "G", "int", context__.to_vec());
            G = int(0);
            vals_i__ = context__.vals_i("G");
            pos__ = 0;
            G = vals_i__[pos__++];
            check_greater_or_equal(function__, "G", G, 1);
            current_statement_begin__ = 46;
            validate_non_negative_index("N", "S", S);
            context__.validate_dims("data initialization", "N", "int", context__.to_vec(S));
            N = std::vector<int>(S, int(0));
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            size_t N_k_0_max__ = S;
            for (size_t k_0__ = 0; k_0__ < N_k_0_max__; ++k_0__) {
                N[k_0__] = vals_i__[pos__++];
            }
            size_t N_i_0_max__ = S;
            for (size_t i_0__ = 0; i_0__ < N_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "N[i_0__]", N[i_0__], 0);
            }
            current_statement_begin__ = 47;
            validate_non_negative_index("T", "S", S);
            context__.validate_dims("data initialization", "T", "double", context__.to_vec(S));
            T = std::vector<double>(S, double(0));
            vals_r__ = context__.vals_r("T");
            pos__ = 0;
            size_t T_k_0_max__ = S;
            for (size_t k_0__ = 0; k_0__ < T_k_0_max__; ++k_0__) {
                T[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 49;
            validate_non_negative_index("t_array", "(G - 1)", (G - 1));
            context__.validate_dims("data initialization", "t_array", "double", context__.to_vec((G - 1)));
            t_array = std::vector<double>((G - 1), double(0));
            vals_r__ = context__.vals_r("t_array");
            pos__ = 0;
            size_t t_array_k_0_max__ = (G - 1);
            for (size_t k_0__ = 0; k_0__ < t_array_k_0_max__; ++k_0__) {
                t_array[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 51;
            context__.validate_dims("data initialization", "prior_K", "double", context__.to_vec());
            prior_K = double(0);
            vals_r__ = context__.vals_r("prior_K");
            pos__ = 0;
            prior_K = vals_r__[pos__++];
            check_greater_or_equal(function__, "prior_K", prior_K, 0);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 55;
            validate_non_negative_index("rho", "G", G);
            num_params_r__ += (1 * G);
            current_statement_begin__ = 57;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_logistic_start_at_1() { }
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
        current_statement_begin__ = 55;
        if (!(context__.contains_r("rho")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable rho missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("rho");
        pos__ = 0U;
        validate_non_negative_index("rho", "G", G);
        context__.validate_dims("parameter initialization", "rho", "double", context__.to_vec(G));
        std::vector<double> rho(G, double(0));
        size_t rho_k_0_max__ = G;
        for (size_t k_0__ = 0; k_0__ < rho_k_0_max__; ++k_0__) {
            rho[k_0__] = vals_r__[pos__++];
        }
        size_t rho_i_0_max__ = G;
        for (size_t i_0__ = 0; i_0__ < rho_i_0_max__; ++i_0__) {
            try {
                writer__.scalar_unconstrain(rho[i_0__]);
            } catch (const std::exception& e) {
                stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable rho: ") + e.what()), current_statement_begin__, prog_reader__());
            }
        }
        current_statement_begin__ = 57;
        if (!(context__.contains_r("K")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable K missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("K");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "K", "double", context__.to_vec());
        double K(0);
        K = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(prior_K, K);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable K: ") + e.what()), current_statement_begin__, prog_reader__());
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
            current_statement_begin__ = 55;
            std::vector<local_scalar_t__> rho;
            size_t rho_d_0_max__ = G;
            rho.reserve(rho_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < rho_d_0_max__; ++d_0__) {
                if (jacobian__)
                    rho.push_back(in__.scalar_constrain(lp__));
                else
                    rho.push_back(in__.scalar_constrain());
            }
            current_statement_begin__ = 57;
            local_scalar_t__ K;
            (void) K;  // dummy to suppress unused var warning
            if (jacobian__)
                K = in__.scalar_lb_constrain(prior_K, lp__);
            else
                K = in__.scalar_lb_constrain(prior_K);
            // model body
            current_statement_begin__ = 62;
            lp_accum__.add(normal_log(K, prior_K, prior_K));
            current_statement_begin__ = 63;
            for (int i = 1; i <= G; ++i) {
                current_statement_begin__ = 64;
                lp_accum__.add(normal_log(get_base1(rho, i, "rho", 1), 0, 1));
            }
            current_statement_begin__ = 67;
            for (int i = 2; i <= S; ++i) {
                current_statement_begin__ = 68;
                lp_accum__.add(poisson_log(get_base1(N, i, "N", 1), mean_t(get_base1(T, i, "T", 1), get_base1(T, 1, "T", 1), get_base1(N, 1, "N", 1), t_array, rho, K, pstream__)));
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
        names__.push_back("rho");
        names__.push_back("K");
        names__.push_back("N_rep");
        names__.push_back("log_lik");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(G);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back((S - 1));
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back((S - 1));
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
        static const char* function__ = "model_logistic_start_at_1_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        std::vector<double> rho;
        size_t rho_d_0_max__ = G;
        rho.reserve(rho_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < rho_d_0_max__; ++d_0__) {
            rho.push_back(in__.scalar_constrain());
        }
        size_t rho_k_0_max__ = G;
        for (size_t k_0__ = 0; k_0__ < rho_k_0_max__; ++k_0__) {
            vars__.push_back(rho[k_0__]);
        }
        double K = in__.scalar_lb_constrain(prior_K);
        vars__.push_back(K);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 73;
            validate_non_negative_index("N_rep", "(S - 1)", (S - 1));
            std::vector<int> N_rep((S - 1), int(0));
            stan::math::fill(N_rep, std::numeric_limits<int>::min());
            current_statement_begin__ = 74;
            validate_non_negative_index("log_lik", "(S - 1)", (S - 1));
            std::vector<double> log_lik((S - 1), double(0));
            stan::math::initialize(log_lik, DUMMY_VAR__);
            stan::math::fill(log_lik, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 76;
            for (int i = 2; i <= S; ++i) {
                current_statement_begin__ = 77;
                stan::model::assign(log_lik, 
                            stan::model::cons_list(stan::model::index_uni((i - 1)), stan::model::nil_index_list()), 
                            poisson_log(get_base1(N, i, "N", 1), mean_t(get_base1(T, i, "T", 1), get_base1(T, 1, "T", 1), get_base1(N, 1, "N", 1), t_array, rho, K, pstream__)), 
                            "assigning variable log_lik");
                current_statement_begin__ = 78;
                stan::model::assign(N_rep, 
                            stan::model::cons_list(stan::model::index_uni((i - 1)), stan::model::nil_index_list()), 
                            poisson_rng(mean_t(get_base1(T, i, "T", 1), get_base1(T, 1, "T", 1), get_base1(N, 1, "N", 1), t_array, rho, K, pstream__), base_rng__), 
                            "assigning variable N_rep");
            }
            // validate, write generated quantities
            current_statement_begin__ = 73;
            size_t N_rep_i_0_max__ = (S - 1);
            for (size_t i_0__ = 0; i_0__ < N_rep_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "N_rep[i_0__]", N_rep[i_0__], 0);
            }
            size_t N_rep_k_0_max__ = (S - 1);
            for (size_t k_0__ = 0; k_0__ < N_rep_k_0_max__; ++k_0__) {
                vars__.push_back(N_rep[k_0__]);
            }
            current_statement_begin__ = 74;
            size_t log_lik_k_0_max__ = (S - 1);
            for (size_t k_0__ = 0; k_0__ < log_lik_k_0_max__; ++k_0__) {
                vars__.push_back(log_lik[k_0__]);
            }
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
        return "model_logistic_start_at_1";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t rho_k_0_max__ = G;
        for (size_t k_0__ = 0; k_0__ < rho_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "rho" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "K";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        size_t N_rep_k_0_max__ = (S - 1);
        for (size_t k_0__ = 0; k_0__ < N_rep_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "N_rep" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t log_lik_k_0_max__ = (S - 1);
        for (size_t k_0__ = 0; k_0__ < log_lik_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t rho_k_0_max__ = G;
        for (size_t k_0__ = 0; k_0__ < rho_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "rho" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "K";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        size_t N_rep_k_0_max__ = (S - 1);
        for (size_t k_0__ = 0; k_0__ < N_rep_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "N_rep" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t log_lik_k_0_max__ = (S - 1);
        for (size_t k_0__ = 0; k_0__ < log_lik_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
}; // model
}  // namespace
typedef model_logistic_start_at_1_namespace::model_logistic_start_at_1 stan_model;
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
