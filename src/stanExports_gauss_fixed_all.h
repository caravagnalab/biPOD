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
namespace model_gauss_fixed_all_namespace {
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
    reader.add_event(0, 0, "start", "model_gauss_fixed_all");
    reader.add_event(59, 57, "end", "model_gauss_fixed_all");
    return reader;
}
template <typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
mean_t(const T0__& n0,
           const T1__& lambda,
           const T2__& mu,
           const T3__& t, std::ostream* pstream__) {
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
        local_scalar_t__ m(DUMMY_VAR__);
        (void) m;  // dummy to suppress unused var warning
        stan::math::initialize(m, DUMMY_VAR__);
        stan::math::fill(m, DUMMY_VAR__);
        current_statement_begin__ = 5;
        stan::math::assign(m, (n0 * stan::math::exp((t * (lambda - mu)))));
        current_statement_begin__ = 6;
        return stan::math::promote_scalar<fun_return_scalar_t__>(m);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct mean_t_functor__ {
    template <typename T0__, typename T1__, typename T2__, typename T3__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
    operator()(const T0__& n0,
           const T1__& lambda,
           const T2__& mu,
           const T3__& t, std::ostream* pstream__) const {
        return mean_t(n0, lambda, mu, t, pstream__);
    }
};
template <typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
sigma_t(const T0__& n0,
            const T1__& lambda,
            const T2__& mu,
            const T3__& t, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 10;
        local_scalar_t__ s(DUMMY_VAR__);
        (void) s;  // dummy to suppress unused var warning
        stan::math::initialize(s, DUMMY_VAR__);
        stan::math::fill(s, DUMMY_VAR__);
        current_statement_begin__ = 11;
        local_scalar_t__ w(DUMMY_VAR__);
        (void) w;  // dummy to suppress unused var warning
        stan::math::initialize(w, DUMMY_VAR__);
        stan::math::fill(w, DUMMY_VAR__);
        stan::math::assign(w,(lambda - mu));
        current_statement_begin__ = 13;
        if (as_bool(logical_eq(lambda, mu))) {
            current_statement_begin__ = 14;
            stan::math::assign(s, (((n0 * 2) * lambda) * t));
        } else {
            current_statement_begin__ = 16;
            stan::math::assign(s, (((n0 * ((lambda + mu) / w)) * stan::math::exp((w * t))) * (stan::math::exp((w * t)) - 1)));
        }
        current_statement_begin__ = 18;
        return stan::math::promote_scalar<fun_return_scalar_t__>(stan::math::sqrt(s));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct sigma_t_functor__ {
    template <typename T0__, typename T1__, typename T2__, typename T3__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
    operator()(const T0__& n0,
            const T1__& lambda,
            const T2__& mu,
            const T3__& t, std::ostream* pstream__) const {
        return sigma_t(n0, lambda, mu, t, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_gauss_fixed_all
  : public stan::model::model_base_crtp<model_gauss_fixed_all> {
private:
        int S;
        double n0;
        std::vector<double> N;
        std::vector<double> T;
public:
    model_gauss_fixed_all(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_gauss_fixed_all(stan::io::var_context& context__,
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
        static const char* function__ = "model_gauss_fixed_all_namespace::model_gauss_fixed_all";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 23;
            context__.validate_dims("data initialization", "S", "int", context__.to_vec());
            S = int(0);
            vals_i__ = context__.vals_i("S");
            pos__ = 0;
            S = vals_i__[pos__++];
            check_greater_or_equal(function__, "S", S, 2);
            current_statement_begin__ = 24;
            context__.validate_dims("data initialization", "n0", "double", context__.to_vec());
            n0 = double(0);
            vals_r__ = context__.vals_r("n0");
            pos__ = 0;
            n0 = vals_r__[pos__++];
            check_greater_or_equal(function__, "n0", n0, 0);
            current_statement_begin__ = 26;
            validate_non_negative_index("N", "S", S);
            context__.validate_dims("data initialization", "N", "double", context__.to_vec(S));
            N = std::vector<double>(S, double(0));
            vals_r__ = context__.vals_r("N");
            pos__ = 0;
            size_t N_k_0_max__ = S;
            for (size_t k_0__ = 0; k_0__ < N_k_0_max__; ++k_0__) {
                N[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 27;
            validate_non_negative_index("T", "S", S);
            context__.validate_dims("data initialization", "T", "double", context__.to_vec(S));
            T = std::vector<double>(S, double(0));
            vals_r__ = context__.vals_r("T");
            pos__ = 0;
            size_t T_k_0_max__ = S;
            for (size_t k_0__ = 0; k_0__ < T_k_0_max__; ++k_0__) {
                T[k_0__] = vals_r__[pos__++];
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 31;
            num_params_r__ += 1;
            current_statement_begin__ = 32;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_gauss_fixed_all() { }
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
        current_statement_begin__ = 31;
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
        current_statement_begin__ = 32;
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
            current_statement_begin__ = 31;
            local_scalar_t__ lambda;
            (void) lambda;  // dummy to suppress unused var warning
            if (jacobian__)
                lambda = in__.scalar_lb_constrain(0, lp__);
            else
                lambda = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 32;
            local_scalar_t__ mu;
            (void) mu;  // dummy to suppress unused var warning
            if (jacobian__)
                mu = in__.scalar_lb_constrain(0, lp__);
            else
                mu = in__.scalar_lb_constrain(0);
            // model body
            current_statement_begin__ = 36;
            lp_accum__.add(inv_gamma_log<propto__>(lambda, 3, 3));
            current_statement_begin__ = 37;
            lp_accum__.add(inv_gamma_log<propto__>(mu, 3, 3));
            current_statement_begin__ = 39;
            for (int i = 1; i <= S; ++i) {
                current_statement_begin__ = 40;
                if (as_bool(logical_eq(i, 1))) {
                    current_statement_begin__ = 41;
                    lp_accum__.add(normal_log<propto__>(get_base1(N, i, "N", 1), mean_t(n0, lambda, mu, get_base1(T, 1, "T", 1), pstream__), sigma_t(n0, lambda, mu, get_base1(T, 1, "T", 1), pstream__)));
                } else {
                    current_statement_begin__ = 43;
                    lp_accum__.add(normal_log<propto__>(get_base1(N, i, "N", 1), mean_t(get_base1(N, (i - 1), "N", 1), lambda, mu, (get_base1(T, i, "T", 1) - get_base1(T, (i - 1), "T", 1)), pstream__), sigma_t(get_base1(N, (i - 1), "N", 1), lambda, mu, (get_base1(T, i, "T", 1) - get_base1(T, (i - 1), "T", 1)), pstream__)));
                }
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
        names__.push_back("N_rep");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(S);
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
        static const char* function__ = "model_gauss_fixed_all_namespace::write_array";
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
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 49;
            validate_non_negative_index("N_rep", "S", S);
            std::vector<double> N_rep(S, double(0));
            stan::math::initialize(N_rep, DUMMY_VAR__);
            stan::math::fill(N_rep, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 50;
            for (int i = 1; i <= S; ++i) {
                current_statement_begin__ = 51;
                if (as_bool(logical_eq(i, 1))) {
                    current_statement_begin__ = 52;
                    stan::model::assign(N_rep, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                normal_rng(mean_t(n0, lambda, mu, get_base1(T, 1, "T", 1), pstream__), sigma_t(n0, lambda, mu, get_base1(T, 1, "T", 1), pstream__), base_rng__), 
                                "assigning variable N_rep");
                } else {
                    current_statement_begin__ = 54;
                    stan::model::assign(N_rep, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                normal_rng(mean_t(get_base1(N, (i - 1), "N", 1), lambda, mu, (get_base1(T, i, "T", 1) - get_base1(T, (i - 1), "T", 1)), pstream__), sigma_t(get_base1(N, (i - 1), "N", 1), lambda, mu, (get_base1(T, i, "T", 1) - get_base1(T, (i - 1), "T", 1)), pstream__), base_rng__), 
                                "assigning variable N_rep");
                }
            }
            // validate, write generated quantities
            current_statement_begin__ = 49;
            size_t N_rep_k_0_max__ = S;
            for (size_t k_0__ = 0; k_0__ < N_rep_k_0_max__; ++k_0__) {
                vars__.push_back(N_rep[k_0__]);
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
        return "model_gauss_fixed_all";
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
        }
        if (!include_gqs__) return;
        size_t N_rep_k_0_max__ = S;
        for (size_t k_0__ = 0; k_0__ < N_rep_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "N_rep" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
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
        }
        if (!include_gqs__) return;
        size_t N_rep_k_0_max__ = S;
        for (size_t k_0__ = 0; k_0__ < N_rep_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "N_rep" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
}; // model
}  // namespace
typedef model_gauss_fixed_all_namespace::model_gauss_fixed_all stan_model;
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
