
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
static constexpr std::array<const char*, 35> locations_array__ = 
{" (found before start of program)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 23, column 4 to column 38)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 28, column 13 to column 24)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 28, column 4 to column 26)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 29, column 4 to column 13)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 31, column 8 to column 40)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 30, column 27 to line 32, column 5)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 30, column 4 to line 32, column 5)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 36, column 8 to column 46)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 35, column 23 to line 37, column 5)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 35, column 4 to line 37, column 5)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 43, column 8 to column 65)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 47, column 12 to column 50)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 48, column 12 to column 92)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 46, column 29 to line 51, column 9)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 46, column 8 to line 51, column 9)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 40, column 26 to line 52, column 5)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 40, column 4 to line 52, column 5)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 2, column 4 to column 27)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 3, column 4 to column 27)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 4, column 4 to column 26)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 5, column 4 to column 35)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 6, column 4 to column 26)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 8, column 25 to column 33)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 8, column 4 to column 35)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 9, column 24 to column 32)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 9, column 4 to column 34)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 10, column 4 to column 23)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 12, column 4 to column 14)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 13, column 4 to column 14)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 14, column 4 to column 22)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 18, column 3 to column 20)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 18, column 21 to column 38)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 19, column 3 to column 30)",
 " (in '/hpc/home/kk319/GitRepo/Proj_CombEffect_STARRseq/notebooks/12_model_sliding/model.stan', line 23, column 29 to column 36)"};



class model_model final : public model_base_crtp<model_model> {

 private:
  int NUM_DATA;
  int NUM_FRAG;
  int NUM_SEG;
  int NUM_SEG_PER_FRAG;
  int NUM_REP;
  std::vector<int> OUTPUT;
  std::vector<int> INPUT;
  double PHI;
  int L_DNA;
  int L_RNA;
  double shrinkageVar;
  double l_rna;
  double l_dna;
  double LIB_RATIO; 
  
 
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
      current_statement__ = 18;
      context__.validate_dims("data initialization","NUM_DATA","int",
           std::vector<size_t>{});
      NUM_DATA = std::numeric_limits<int>::min();
      
      current_statement__ = 18;
      NUM_DATA = context__.vals_i("NUM_DATA")[(1 - 1)];
      current_statement__ = 18;
      check_greater_or_equal(function__, "NUM_DATA", NUM_DATA, 1);
      current_statement__ = 19;
      context__.validate_dims("data initialization","NUM_FRAG","int",
           std::vector<size_t>{});
      NUM_FRAG = std::numeric_limits<int>::min();
      
      current_statement__ = 19;
      NUM_FRAG = context__.vals_i("NUM_FRAG")[(1 - 1)];
      current_statement__ = 19;
      check_greater_or_equal(function__, "NUM_FRAG", NUM_FRAG, 1);
      current_statement__ = 20;
      context__.validate_dims("data initialization","NUM_SEG","int",
           std::vector<size_t>{});
      NUM_SEG = std::numeric_limits<int>::min();
      
      current_statement__ = 20;
      NUM_SEG = context__.vals_i("NUM_SEG")[(1 - 1)];
      current_statement__ = 20;
      check_greater_or_equal(function__, "NUM_SEG", NUM_SEG, 1);
      current_statement__ = 21;
      context__.validate_dims("data initialization","NUM_SEG_PER_FRAG","int",
           std::vector<size_t>{});
      NUM_SEG_PER_FRAG = std::numeric_limits<int>::min();
      
      current_statement__ = 21;
      NUM_SEG_PER_FRAG = context__.vals_i("NUM_SEG_PER_FRAG")[(1 - 1)];
      current_statement__ = 21;
      check_greater_or_equal(function__, "NUM_SEG_PER_FRAG",
                             NUM_SEG_PER_FRAG, 1);
      current_statement__ = 22;
      context__.validate_dims("data initialization","NUM_REP","int",
           std::vector<size_t>{});
      NUM_REP = std::numeric_limits<int>::min();
      
      current_statement__ = 22;
      NUM_REP = context__.vals_i("NUM_REP")[(1 - 1)];
      current_statement__ = 22;
      check_greater_or_equal(function__, "NUM_REP", NUM_REP, 1);
      current_statement__ = 23;
      validate_non_negative_index("OUTPUT", "NUM_DATA", NUM_DATA);
      current_statement__ = 24;
      context__.validate_dims("data initialization","OUTPUT","int",
           std::vector<size_t>{static_cast<size_t>(NUM_DATA)});
      OUTPUT = std::vector<int>(NUM_DATA, std::numeric_limits<int>::min());
      
      current_statement__ = 24;
      OUTPUT = context__.vals_i("OUTPUT");
      current_statement__ = 24;
      for (int sym1__ = 1; sym1__ <= NUM_DATA; ++sym1__) {
        current_statement__ = 24;
        check_greater_or_equal(function__, "OUTPUT[sym1__]",
                               OUTPUT[(sym1__ - 1)], 0);
      }
      current_statement__ = 25;
      validate_non_negative_index("INPUT", "NUM_DATA", NUM_DATA);
      current_statement__ = 26;
      context__.validate_dims("data initialization","INPUT","int",
           std::vector<size_t>{static_cast<size_t>(NUM_DATA)});
      INPUT = std::vector<int>(NUM_DATA, std::numeric_limits<int>::min());
      
      current_statement__ = 26;
      INPUT = context__.vals_i("INPUT");
      current_statement__ = 26;
      for (int sym1__ = 1; sym1__ <= NUM_DATA; ++sym1__) {
        current_statement__ = 26;
        check_greater_or_equal(function__, "INPUT[sym1__]",
                               INPUT[(sym1__ - 1)], 1);
      }
      current_statement__ = 27;
      context__.validate_dims("data initialization","PHI","double",
           std::vector<size_t>{});
      PHI = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 27;
      PHI = context__.vals_r("PHI")[(1 - 1)];
      current_statement__ = 27;
      check_greater_or_equal(function__, "PHI", PHI, 0);
      current_statement__ = 28;
      context__.validate_dims("data initialization","L_DNA","int",
           std::vector<size_t>{});
      L_DNA = std::numeric_limits<int>::min();
      
      current_statement__ = 28;
      L_DNA = context__.vals_i("L_DNA")[(1 - 1)];
      current_statement__ = 29;
      context__.validate_dims("data initialization","L_RNA","int",
           std::vector<size_t>{});
      L_RNA = std::numeric_limits<int>::min();
      
      current_statement__ = 29;
      L_RNA = context__.vals_i("L_RNA")[(1 - 1)];
      current_statement__ = 30;
      context__.validate_dims("data initialization","shrinkageVar","double",
           std::vector<size_t>{});
      shrinkageVar = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 30;
      shrinkageVar = context__.vals_r("shrinkageVar")[(1 - 1)];
      current_statement__ = 31;
      l_rna = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 31;
      l_rna = L_RNA;
      current_statement__ = 32;
      l_dna = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 32;
      l_dna = L_DNA;
      current_statement__ = 33;
      LIB_RATIO = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 33;
      LIB_RATIO = (l_rna / l_dna);
      current_statement__ = 34;
      validate_non_negative_index("gamma", "NUM_SEG", NUM_SEG);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = NUM_SEG;
    
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
      std::vector<local_scalar_t__> gamma;
      gamma = std::vector<local_scalar_t__>(NUM_SEG, DUMMY_VAR__);
      
      current_statement__ = 1;
      gamma = in__.template read_constrain_lb<std::vector<local_scalar_t__>, jacobian__>(
                0.0001, lp__, NUM_SEG);
      {
        current_statement__ = 2;
        validate_non_negative_index("PSA", "NUM_SEG + 1", (NUM_SEG + 1));
        std::vector<local_scalar_t__> PSA;
        PSA = std::vector<local_scalar_t__>((NUM_SEG + 1), DUMMY_VAR__);
        
        current_statement__ = 4;
        assign(PSA, 0, "assigning variable PSA", index_uni(1));
        current_statement__ = 7;
        for (int i = 2; i <= (NUM_SEG + 1); ++i) {
          current_statement__ = 5;
          assign(PSA,
            (rvalue(PSA, "PSA", index_uni((i - 1))) +
              stan::math::log(rvalue(gamma, "gamma", index_uni((i - 1))))),
            "assigning variable PSA", index_uni(i));
        }
        current_statement__ = 10;
        for (int i = 1; i <= NUM_SEG; ++i) {
          current_statement__ = 8;
          lp_accum__.add(
            lognormal_lpdf<propto__>(rvalue(gamma, "gamma", index_uni(i)), 0,
              shrinkageVar));
        }
        current_statement__ = 17;
        for (int i = 1; i <= NUM_FRAG; ++i) {
          local_scalar_t__ gammaProd;
          gammaProd = DUMMY_VAR__;
          
          current_statement__ = 11;
          gammaProd = stan::math::exp(
                        (rvalue(PSA, "PSA",
                           index_uni((i + NUM_SEG_PER_FRAG))) -
                          rvalue(PSA, "PSA", index_uni(i))));
          current_statement__ = 15;
          for (int k = 1; k <= NUM_REP; ++k) {
            int idx_data;
            idx_data = std::numeric_limits<int>::min();
            
            current_statement__ = 12;
            idx_data = ((NUM_FRAG * (i - 1)) + k);
            current_statement__ = 13;
            lp_accum__.add(
              neg_binomial_2_lpmf<propto__>(
                rvalue(OUTPUT, "OUTPUT", index_uni(idx_data)),
                ((rvalue(INPUT, "INPUT", index_uni(idx_data)) * LIB_RATIO) *
                  gammaProd), PHI));
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
      std::vector<double> gamma;
      gamma = std::vector<double>(NUM_SEG, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 1;
      gamma = in__.template read_constrain_lb<std::vector<local_scalar_t__>, jacobian__>(
                0.0001, lp__, NUM_SEG);
      for (int sym1__ = 1; sym1__ <= NUM_SEG; ++sym1__) {
        vars__.emplace_back(gamma[(sym1__ - 1)]);
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
      std::vector<double> gamma;
      gamma = std::vector<double>(NUM_SEG, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 1;
      gamma = context__.vals_r("gamma");
      std::vector<double> gamma_free__;
      gamma_free__ = std::vector<double>(NUM_SEG, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 1;
      for (int sym1__ = 1; sym1__ <= NUM_SEG; ++sym1__) {
        current_statement__ = 1;
        assign(gamma_free__,
          stan::math::lb_free(gamma[(sym1__ - 1)], 0.0001),
          "assigning variable gamma_free__", index_uni(sym1__));
      }
      for (int sym1__ = 1; sym1__ <= NUM_SEG; ++sym1__) {
        vars__.emplace_back(gamma_free__[(sym1__ - 1)]);
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__ = std::vector<std::string>{"gamma"};
    
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{
                                                                   static_cast<size_t>(NUM_SEG)
                                                                   }};
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    for (int sym1__ = 1; sym1__ <= NUM_SEG; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "gamma" + '.' + std::to_string(sym1__));
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
    
    for (int sym1__ = 1; sym1__ <= NUM_SEG; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "gamma" + '.' + std::to_string(sym1__));
      } 
    }
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"gamma\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(NUM_SEG) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"}]");
    
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"gamma\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(NUM_SEG) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"}]");
    
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


