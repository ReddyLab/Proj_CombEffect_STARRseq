
// Code generated by stanc 888c0f6a
#include <stan/model/model_header.hpp>
namespace model_model_namespace {

using stan::io::dump;
using stan::model::assign;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using namespace stan::math;


stan::math::profile_map profiles__;
static constexpr std::array<const char*, 27> locations_array__ = 
{" (found before start of program)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 14, column 4 to column 29)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 25, column 8 to column 35)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 24, column 23 to line 26, column 5)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 24, column 4 to line 26, column 5)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 32, column 8 to column 28)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 34, column 12 to column 33)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 33, column 46 to line 35, column 9)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 33, column 8 to line 35, column 9)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 38, column 8 to column 21)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 40, column 12 to column 44)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 41, column 12 to column 81)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 39, column 27 to line 44, column 9)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 39, column 8 to line 44, column 9)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 29, column 24 to line 45, column 5)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 29, column 4 to line 45, column 5)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 2, column 4 to column 25)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 3, column 4 to column 25)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 4, column 4 to column 24)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 5, column 4 to column 33)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 6, column 4 to column 24)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 8, column 25 to column 31)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 8, column 4 to column 33)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 9, column 24 to column 30)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 9, column 4 to column 32)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 10, column 4 to column 23)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/model_kmer_sim/model.stan', line 14, column 22 to column 27)"};



class model_model final : public model_base_crtp<model_model> {

 private:
  int N_DATA;
  int N_FRAG;
  int N_SEG;
  int N_SEG_PER_FRAG;
  int N_REP;
  std::vector<int> OUTPUT;
  std::vector<int> INPUT;
  double PHI; 
  
 
 public:
  ~model_model() { }
  
  inline std::string model_name() const final { return "model_model"; }

  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 888c0f6a", "stancflags = "};
  }
  
  
  model_model(stan::io::var_context& context__,
              unsigned int random_seed__ = 0,
              std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static constexpr const char* function__ = "model_model_namespace::model_model";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 16;
      context__.validate_dims("data initialization","N_DATA","int",
           std::vector<size_t>{});
      N_DATA = std::numeric_limits<int>::min();
      
      current_statement__ = 16;
      N_DATA = context__.vals_i("N_DATA")[(1 - 1)];
      current_statement__ = 16;
      check_greater_or_equal(function__, "N_DATA", N_DATA, 1);
      current_statement__ = 17;
      context__.validate_dims("data initialization","N_FRAG","int",
           std::vector<size_t>{});
      N_FRAG = std::numeric_limits<int>::min();
      
      current_statement__ = 17;
      N_FRAG = context__.vals_i("N_FRAG")[(1 - 1)];
      current_statement__ = 17;
      check_greater_or_equal(function__, "N_FRAG", N_FRAG, 1);
      current_statement__ = 18;
      context__.validate_dims("data initialization","N_SEG","int",
           std::vector<size_t>{});
      N_SEG = std::numeric_limits<int>::min();
      
      current_statement__ = 18;
      N_SEG = context__.vals_i("N_SEG")[(1 - 1)];
      current_statement__ = 18;
      check_greater_or_equal(function__, "N_SEG", N_SEG, 1);
      current_statement__ = 19;
      context__.validate_dims("data initialization","N_SEG_PER_FRAG","int",
           std::vector<size_t>{});
      N_SEG_PER_FRAG = std::numeric_limits<int>::min();
      
      current_statement__ = 19;
      N_SEG_PER_FRAG = context__.vals_i("N_SEG_PER_FRAG")[(1 - 1)];
      current_statement__ = 19;
      check_greater_or_equal(function__, "N_SEG_PER_FRAG", N_SEG_PER_FRAG, 1);
      current_statement__ = 20;
      context__.validate_dims("data initialization","N_REP","int",
           std::vector<size_t>{});
      N_REP = std::numeric_limits<int>::min();
      
      current_statement__ = 20;
      N_REP = context__.vals_i("N_REP")[(1 - 1)];
      current_statement__ = 20;
      check_greater_or_equal(function__, "N_REP", N_REP, 1);
      current_statement__ = 21;
      validate_non_negative_index("OUTPUT", "N_DATA", N_DATA);
      current_statement__ = 22;
      context__.validate_dims("data initialization","OUTPUT","int",
           std::vector<size_t>{static_cast<size_t>(N_DATA)});
      OUTPUT = std::vector<int>(N_DATA, std::numeric_limits<int>::min());
      
      current_statement__ = 22;
      OUTPUT = context__.vals_i("OUTPUT");
      current_statement__ = 22;
      for (int sym1__ = 1; sym1__ <= N_DATA; ++sym1__) {
        current_statement__ = 22;
        check_greater_or_equal(function__, "OUTPUT[sym1__]",
                               OUTPUT[(sym1__ - 1)], 0);
      }
      current_statement__ = 23;
      validate_non_negative_index("INPUT", "N_DATA", N_DATA);
      current_statement__ = 24;
      context__.validate_dims("data initialization","INPUT","int",
           std::vector<size_t>{static_cast<size_t>(N_DATA)});
      INPUT = std::vector<int>(N_DATA, std::numeric_limits<int>::min());
      
      current_statement__ = 24;
      INPUT = context__.vals_i("INPUT");
      current_statement__ = 24;
      for (int sym1__ = 1; sym1__ <= N_DATA; ++sym1__) {
        current_statement__ = 24;
        check_greater_or_equal(function__, "INPUT[sym1__]",
                               INPUT[(sym1__ - 1)], 1);
      }
      current_statement__ = 25;
      context__.validate_dims("data initialization","PHI","double",
           std::vector<size_t>{});
      PHI = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 25;
      PHI = context__.vals_r("PHI")[(1 - 1)];
      current_statement__ = 25;
      check_greater_or_equal(function__, "PHI", PHI, 0);
      current_statement__ = 26;
      validate_non_negative_index("SEG", "N_SEG", N_SEG);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = N_SEG;
    
  }
  
  template <bool propto__, bool jacobian__ , typename VecR, typename VecI, 
  stan::require_vector_like_t<VecR>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr> 
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "model_model_namespace::log_prob";
    (void) function__;  // suppress unused var warning
    
    try {
      std::vector<local_scalar_t__> SEG;
      SEG = std::vector<local_scalar_t__>(N_SEG, DUMMY_VAR__);
      
      current_statement__ = 1;
      SEG = in__.template read_constrain_lb<std::vector<local_scalar_t__>, jacobian__>(
              0, lp__, N_SEG);
      {
        current_statement__ = 4;
        for (int j = 1; j <= N_SEG; ++j) {
          current_statement__ = 2;
          lp_accum__.add(
            lognormal_lpdf<propto__>(rvalue(SEG, "SEG", index_uni(j)), 0,
              0.3));
        }
        current_statement__ = 15;
        for (int i = 1; i <= N_FRAG; ++i) {
          local_scalar_t__ seg_effect;
          seg_effect = DUMMY_VAR__;
          
          current_statement__ = 5;
          seg_effect = 1;
          current_statement__ = 8;
          for (int j = i; j <= ((i + N_SEG_PER_FRAG) - 1); ++j) {
            current_statement__ = 6;
            seg_effect = (seg_effect * rvalue(SEG, "SEG", index_uni(j)));
          }
          int idx_data;
          idx_data = std::numeric_limits<int>::min();
          
          current_statement__ = 13;
          for (int k = 1; k <= N_REP; ++k) {
            current_statement__ = 10;
            idx_data = ((N_FRAG * (i - 1)) + k);
            current_statement__ = 11;
            lp_accum__.add(
              neg_binomial_2_lpmf<propto__>(
                rvalue(OUTPUT, "OUTPUT", index_uni(idx_data)),
                (rvalue(INPUT, "INPUT", index_uni(idx_data)) * seg_effect),
                PHI));
          }
        }
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, 
  stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, 
  stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr> 
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.resize(0);
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    static constexpr bool propto__ = true;
    (void) propto__;
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    int current_statement__ = 0; 
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    constexpr bool jacobian__ = false;
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "model_model_namespace::write_array";
    (void) function__;  // suppress unused var warning
    
    try {
      std::vector<double> SEG;
      SEG = std::vector<double>(N_SEG, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 1;
      SEG = in__.template read_constrain_lb<std::vector<local_scalar_t__>, jacobian__>(
              0, lp__, N_SEG);
      for (int sym1__ = 1; sym1__ <= N_SEG; ++sym1__) {
        vars__.emplace_back(SEG[(sym1__ - 1)]);
      }
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, 
  stan::require_std_vector_t<VecVar>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr> 
  inline void transform_inits_impl(const stan::io::var_context& context__,
                                   VecI& params_i__, VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.clear();
    vars__.reserve(num_params_r__);
    int current_statement__ = 0; 
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      std::vector<double> SEG;
      SEG = std::vector<double>(N_SEG, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 1;
      SEG = context__.vals_r("SEG");
      std::vector<double> SEG_free__;
      SEG_free__ = std::vector<double>(N_SEG, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 1;
      for (int sym1__ = 1; sym1__ <= N_SEG; ++sym1__) {
        current_statement__ = 1;
        assign(SEG_free__, stan::math::lb_free(SEG[(sym1__ - 1)], 0),
          "assigning variable SEG_free__", index_uni(sym1__));
      }
      for (int sym1__ = 1; sym1__ <= N_SEG; ++sym1__) {
        vars__.emplace_back(SEG_free__[(sym1__ - 1)]);
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__ = std::vector<std::string>{"SEG"};
    
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{
                                                                   static_cast<size_t>(N_SEG)
                                                                   }};
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    for (int sym1__ = 1; sym1__ <= N_SEG; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "SEG" + '.' + std::to_string(sym1__));
      } 
    }
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    for (int sym1__ = 1; sym1__ <= N_SEG; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "SEG" + '.' + std::to_string(sym1__));
      } 
    }
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"SEG\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(N_SEG) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"}]");
    
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"SEG\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(N_SEG) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"}]");
    
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      std::vector<double> vars_vec;
      vars_vec.reserve(vars.size());
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>>(
        vars_vec.data(), vars_vec.size());
    }

    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      write_array_impl(base_rng, params_r, params_i, vars,
       emit_transformed_parameters, emit_generated_quantities, pstream);
    }

    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }

    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }


    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec;
      params_r_vec.reserve(params_r.size());
      std::vector<int> params_i;
      transform_inits_impl(context, params_i, params_r_vec, pstream);
      params_r = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>>(
        params_r_vec.data(), params_r_vec.size());
    }
    inline void transform_inits(const stan::io::var_context& context,
                                std::vector<int>& params_i,
                                std::vector<double>& vars,
                                std::ostream* pstream = nullptr) const final {
      transform_inits_impl(context, params_i, vars, pstream);
    }

};
}
using stan_model = model_model_namespace::model_model;

#ifndef USING_R

// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

stan::math::profile_map& get_stan_profile_data() {
  return model_model_namespace::profiles__;
}

#endif


