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
namespace model_exponential_with_changing_points_namespace {
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
    reader.add_event(0, 0, "start", "model_exponential_with_changing_points");
    reader.add_event(76, 74, "end", "model_exponential_with_changing_points");
    return reader;
}
template <typename T1__, typename T2__, typename T3__, typename T4__>
typename boost::math::tools::promote_args<T1__, T2__, T3__, T4__>::type
mean_t(const int& n0,
           const T1__& t0,
           const T2__& t,
           const std::vector<T3__>& changing_times,
           const std::vector<T4__>& rho_array, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T1__, T2__, T3__, T4__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 5;
        local_scalar_t__ res(DUMMY_VAR__);
        (void) res;  // dummy to suppress unused var warning
        stan::math::initialize(res, DUMMY_VAR__);
        stan::math::fill(res, DUMMY_VAR__);
        stan::math::assign(res,n0);
        current_statement_begin__ = 6;
        int J(0);
        (void) J;  // dummy to suppress unused var warning
        stan::math::fill(J, std::numeric_limits<int>::min());
        stan::math::assign(J,num_elements(changing_times));
        current_statement_begin__ = 8;
        if (as_bool(logical_lte(t, get_base1(changing_times, 1, "changing_times", 1)))) {
            current_statement_begin__ = 9;
            return stan::math::promote_scalar<fun_return_scalar_t__>((n0 * stan::math::exp((get_base1(rho_array, 1, "rho_array", 1) * (t - t0)))));
        } else {
            current_statement_begin__ = 11;
            stan::math::assign(res, (res * stan::math::exp((get_base1(rho_array, 1, "rho_array", 1) * (get_base1(changing_times, 1, "changing_times", 1) - t0)))));
            current_statement_begin__ = 12;
            for (int j = 2; j <= J; ++j) {
                current_statement_begin__ = 13;
                if (as_bool(logical_gte(t, get_base1(changing_times, j, "changing_times", 1)))) {
                    current_statement_begin__ = 14;
                    stan::math::assign(res, (res * stan::math::exp((get_base1(rho_array, j, "rho_array", 1) * (get_base1(changing_times, j, "changing_times", 1) - get_base1(changing_times, (j - 1), "changing_times", 1))))));
                } else {
                    current_statement_begin__ = 16;
                    stan::math::assign(res, (res * stan::math::exp((get_base1(rho_array, j, "rho_array", 1) * (t - get_base1(changing_times, (j - 1), "changing_times", 1))))));
                    current_statement_begin__ = 17;
                    return stan::math::promote_scalar<fun_return_scalar_t__>(res);
                }
            }
            current_statement_begin__ = 20;
            stan::math::assign(res, (res * stan::math::exp((get_base1(rho_array, (J + 1), "rho_array", 1) * (t - get_base1(changing_times, J, "changing_times", 1))))));
            current_statement_begin__ = 21;
            return stan::math::promote_scalar<fun_return_scalar_t__>(res);
        }
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct mean_t_functor__ {
    template <typename T1__, typename T2__, typename T3__, typename T4__>
        typename boost::math::tools::promote_args<T1__, T2__, T3__, T4__>::type
    operator()(const int& n0,
           const T1__& t0,
           const T2__& t,
           const std::vector<T3__>& changing_times,
           const std::vector<T4__>& rho_array, std::ostream* pstream__) const {
        return mean_t(n0, t0, t, changing_times, rho_array, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_exponential_with_changing_points
  : public stan::model::model_base_crtp<model_exponential_with_changing_points> {
private:
        int S;
        int G;
        std::vector<int> N;
        std::vector<double> T;
        std::vector<double> changing_times_prior;
        double dt;
public:
    model_exponential_with_changing_points(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_exponential_with_changing_points(stan::io::var_context& context__,
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
        static const char* function__ = "model_exponential_with_changing_points_namespace::model_exponential_with_changing_points";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 28;
            context__.validate_dims("data initialization", "S", "int", context__.to_vec());
            S = int(0);
            vals_i__ = context__.vals_i("S");
            pos__ = 0;
            S = vals_i__[pos__++];
            check_greater_or_equal(function__, "S", S, 1);
            current_statement_begin__ = 29;
            context__.validate_dims("data initialization", "G", "int", context__.to_vec());
            G = int(0);
            vals_i__ = context__.vals_i("G");
            pos__ = 0;
            G = vals_i__[pos__++];
            check_greater_or_equal(function__, "G", G, 1);
            current_statement_begin__ = 31;
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
            current_statement_begin__ = 32;
            validate_non_negative_index("T", "S", S);
            context__.validate_dims("data initialization", "T", "double", context__.to_vec(S));
            T = std::vector<double>(S, double(0));
            vals_r__ = context__.vals_r("T");
            pos__ = 0;
            size_t T_k_0_max__ = S;
            for (size_t k_0__ = 0; k_0__ < T_k_0_max__; ++k_0__) {
                T[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 34;
            validate_non_negative_index("changing_times_prior", "(G - 1)", (G - 1));
            context__.validate_dims("data initialization", "changing_times_prior", "double", context__.to_vec((G - 1)));
            changing_times_prior = std::vector<double>((G - 1), double(0));
            vals_r__ = context__.vals_r("changing_times_prior");
            pos__ = 0;
            size_t changing_times_prior_k_0_max__ = (G - 1);
            for (size_t k_0__ = 0; k_0__ < changing_times_prior_k_0_max__; ++k_0__) {
                changing_times_prior[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 35;
            context__.validate_dims("data initialization", "dt", "double", context__.to_vec());
            dt = double(0);
            vals_r__ = context__.vals_r("dt");
            pos__ = 0;
            dt = vals_r__[pos__++];
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 39;
            validate_non_negative_index("rho", "G", G);
            num_params_r__ += (1 * G);
            current_statement_begin__ = 40;
            validate_non_negative_index("changing_times_unit", "(G - 1)", (G - 1));
            num_params_r__ += (1 * (G - 1));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_exponential_with_changing_points() { }
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
        current_statement_begin__ = 39;
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
        current_statement_begin__ = 40;
        if (!(context__.contains_r("changing_times_unit")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable changing_times_unit missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("changing_times_unit");
        pos__ = 0U;
        validate_non_negative_index("changing_times_unit", "(G - 1)", (G - 1));
        context__.validate_dims("parameter initialization", "changing_times_unit", "double", context__.to_vec((G - 1)));
        std::vector<double> changing_times_unit((G - 1), double(0));
        size_t changing_times_unit_k_0_max__ = (G - 1);
        for (size_t k_0__ = 0; k_0__ < changing_times_unit_k_0_max__; ++k_0__) {
            changing_times_unit[k_0__] = vals_r__[pos__++];
        }
        size_t changing_times_unit_i_0_max__ = (G - 1);
        for (size_t i_0__ = 0; i_0__ < changing_times_unit_i_0_max__; ++i_0__) {
            try {
                writer__.scalar_lub_unconstrain(-(dt), dt, changing_times_unit[i_0__]);
            } catch (const std::exception& e) {
                stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable changing_times_unit: ") + e.what()), current_statement_begin__, prog_reader__());
            }
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
            current_statement_begin__ = 39;
            std::vector<local_scalar_t__> rho;
            size_t rho_d_0_max__ = G;
            rho.reserve(rho_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < rho_d_0_max__; ++d_0__) {
                if (jacobian__)
                    rho.push_back(in__.scalar_constrain(lp__));
                else
                    rho.push_back(in__.scalar_constrain());
            }
            current_statement_begin__ = 40;
            std::vector<local_scalar_t__> changing_times_unit;
            size_t changing_times_unit_d_0_max__ = (G - 1);
            changing_times_unit.reserve(changing_times_unit_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < changing_times_unit_d_0_max__; ++d_0__) {
                if (jacobian__)
                    changing_times_unit.push_back(in__.scalar_lub_constrain(-(dt), dt, lp__));
                else
                    changing_times_unit.push_back(in__.scalar_lub_constrain(-(dt), dt));
            }
            // transformed parameters
            current_statement_begin__ = 44;
            validate_non_negative_index("changing_times", "(G - 1)", (G - 1));
            std::vector<local_scalar_t__> changing_times((G - 1), local_scalar_t__(0));
            stan::math::initialize(changing_times, DUMMY_VAR__);
            stan::math::fill(changing_times, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 46;
            for (int i = 2; i <= G; ++i) {
                current_statement_begin__ = 47;
                stan::model::assign(changing_times, 
                            stan::model::cons_list(stan::model::index_uni((i - 1)), stan::model::nil_index_list()), 
                            (get_base1(changing_times_prior, (i - 1), "changing_times_prior", 1) + get_base1(changing_times_unit, (i - 1), "changing_times_unit", 1)), 
                            "assigning variable changing_times");
            }
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 44;
            size_t changing_times_k_0_max__ = (G - 1);
            for (size_t k_0__ = 0; k_0__ < changing_times_k_0_max__; ++k_0__) {
                if (stan::math::is_uninitialized(changing_times[k_0__])) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: changing_times" << "[" << k_0__ << "]";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable changing_times: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            current_statement_begin__ = 52;
            for (int i = 1; i <= G; ++i) {
                current_statement_begin__ = 53;
                lp_accum__.add(normal_log(get_base1(rho, i, "rho", 1), 0, 1));
            }
            current_statement_begin__ = 56;
            for (int i = 2; i <= G; ++i) {
                current_statement_begin__ = 57;
                lp_accum__.add(uniform_log(get_base1(changing_times_unit, (i - 1), "changing_times_unit", 1), -(dt), dt));
            }
            current_statement_begin__ = 61;
            for (int i = 1; i <= S; ++i) {
                current_statement_begin__ = 62;
                lp_accum__.add(poisson_log(get_base1(N, i, "N", 1), mean_t(get_base1(N, 1, "N", 1), get_base1(T, 1, "T", 1), get_base1(T, i, "T", 1), changing_times, rho, pstream__)));
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
        names__.push_back("changing_times_unit");
        names__.push_back("changing_times");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(G);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back((G - 1));
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back((G - 1));
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
        static const char* function__ = "model_exponential_with_changing_points_namespace::write_array";
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
        std::vector<double> changing_times_unit;
        size_t changing_times_unit_d_0_max__ = (G - 1);
        changing_times_unit.reserve(changing_times_unit_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < changing_times_unit_d_0_max__; ++d_0__) {
            changing_times_unit.push_back(in__.scalar_lub_constrain(-(dt), dt));
        }
        size_t changing_times_unit_k_0_max__ = (G - 1);
        for (size_t k_0__ = 0; k_0__ < changing_times_unit_k_0_max__; ++k_0__) {
            vars__.push_back(changing_times_unit[k_0__]);
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 44;
            validate_non_negative_index("changing_times", "(G - 1)", (G - 1));
            std::vector<double> changing_times((G - 1), double(0));
            stan::math::initialize(changing_times, DUMMY_VAR__);
            stan::math::fill(changing_times, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 46;
            for (int i = 2; i <= G; ++i) {
                current_statement_begin__ = 47;
                stan::model::assign(changing_times, 
                            stan::model::cons_list(stan::model::index_uni((i - 1)), stan::model::nil_index_list()), 
                            (get_base1(changing_times_prior, (i - 1), "changing_times_prior", 1) + get_base1(changing_times_unit, (i - 1), "changing_times_unit", 1)), 
                            "assigning variable changing_times");
            }
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t changing_times_k_0_max__ = (G - 1);
                for (size_t k_0__ = 0; k_0__ < changing_times_k_0_max__; ++k_0__) {
                    vars__.push_back(changing_times[k_0__]);
                }
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
        return "model_exponential_with_changing_points";
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
        size_t changing_times_unit_k_0_max__ = (G - 1);
        for (size_t k_0__ = 0; k_0__ < changing_times_unit_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "changing_times_unit" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t changing_times_k_0_max__ = (G - 1);
            for (size_t k_0__ = 0; k_0__ < changing_times_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "changing_times" << '.' << k_0__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
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
        size_t changing_times_unit_k_0_max__ = (G - 1);
        for (size_t k_0__ = 0; k_0__ < changing_times_unit_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "changing_times_unit" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t changing_times_k_0_max__ = (G - 1);
            for (size_t k_0__ = 0; k_0__ < changing_times_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "changing_times" << '.' << k_0__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_exponential_with_changing_points_namespace::model_exponential_with_changing_points stan_model;
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
