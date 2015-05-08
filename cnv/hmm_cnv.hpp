#define INTERCEPT .32

class hmm_cnv_settings_t{
};

struct sample_info_t{
  string filename;
  string annot;
  float global_alpha;
};

//struct state_t;

class HMM_CNV:public Analyzer{
public:
  static const int INTENSITY_LRR=0;
  static const int INTENSITY_READS=1;
  static const int ANALYSIS_TUMOR=0;
  static const int ANALYSIS_SEQ=1;

  void init(const ptree & pt);
  HMM_CNV(const char * analysis_type);
  ~HMM_CNV();
  void run();
  static const float SMALL=1e-10;
  static const float LARGE=1e10;
  float LOG_SMALL;
  float LOG_LARGE;
  static const float HET_RATIO=.3;
  static const int TUMOR_CN_ELEMENTS = 5;
  static const float SLOPE = (1-INTERCEPT)/2;
private:
  int analysis_type;
#ifdef USE_GPU
  opencl_hmm_cnv_t opencl_info;
  cl_int err;
#endif
  // These are input settings
  string output_dir;
  string pfb_filename;
  string kernel_path;
  string input_filelist;
  // sequencing tuning parameters
  float coverage;
  float coverage_offset;
  float shape;
  float dispersion;
  // tumor tuning paramters
  float * sd_rr_alpha;
  float * sd_rr_loh;
  float * sd_baf_percent;
  float sd_rr_normal;
  float sd_baf_homo_del;
  float global_alpha;
  float normal_cdf;
  int normal_state;
  int homo_deletion_state;

  int intensity_type;
  bool verbose_output;
  MathUtils * math;
  float * penalty_mat; // applies further penalties to transition matrix
  //vector<string> filelist_vec; // list of filenames for input
  vector<sample_info_t> sample_info_vec;
  //vector<float> global_alphas;
  vector<int> alpha_percents;
  string chr; // chr to operate on
  int * pos_vec; // positions in base pairs
  int * gap_vec; // inter-SNP distance in base pairs
  //vector<state_t> state_object_vector; // vector of state objects
  string * rslist; // list of SNP IDs
  int samples; // total samples
  int total_markers;
  int current_sample;
  int observations; // total markers
  int states; // total states
  int iteration;  // current iteration
  float sd_alpha;
  int hmm_matrix_size;
  int trans_matrix_size;
 

  // observed data here
  int * segment_length_vec; // relevant for penncnv-seq
  float *intensity_vec;
  float *baf_vec;
  map<int,float>  pfb_map;
  // these are large matrices MARKERS  * SUBJECTS
  float *forward_mat,*backward_mat, *prob_mat, * prev_prob_mat;
  float * trans_mat; // transition matrix
  //float * trans_mat; // previous transition matrix
  float * emission_cache_mat;
  // these track the number of rescalings needed to prevent underflow
  // dimension is number of markers
  int * forward_rescalings;
  int * backward_rescalings;
  // viterbi variables
  int * greedy_rescalings;
  int * backtrace_vec;
  float * greedy_mat;
  int * bestpath_mat;

  float likelihood_log;
  // needed by GPU
  int * state_tumorcn_vec; // for use on GPU. copy numbers
  bool * state_het_ident;
  bool * state_alpha_ident;
  bool * state_tumor_ident;
  float * state_het_ratio;
  float * state_total_alpha;
  float * state_mu_rr;
  float * state_sd_rr;
  //float * state_baf_het_mean;
  float * state_alpha_prior;
  
  // functions below
  float get_baf_mean(int bac,int cn, float alpha);
  float get_baf_sd(int state,int bac);
  bool load_input_filelist(const string & filelist);
  void load_input_file();
  void init_emission_cache_mat_orig();
  void compute_emission_cache_mat();
  void print_summary();
  void init_state_objects();
  void allocate_global_matrices();
  void allocate_sample_matrices();
  void initialize_global_matrices();
  void initialize_state_parameters();
  void initialize_state_parameter(int id,int tumor_cn,  int percent,  float het_ratio,bool alpha_identifiable, bool het_identifiable, bool tumor_identifiable, float * sd_rr_alpha, float * sd_rr_loh, float  sd_rr_normal);
  void free_sample_matrices();
  void adjust_aneuploidy();
  void forwardbackward();
  void viterbi();
  void viterbi_output();
  void baumwelch();
  float get_rr_offset();
  float get_mean_depth();
  void normalize_across_alphas(float * emission_vec);
  // utilities
  float get_binom_sum(int n, float pop_b, float baf);
  void parse_pfb_file();
  double addMultipleLogs(double * logs, int len);
  void normalize(int len, float * probs, float * normalized);
  //void set_global_alpha();
  //ofstream ofs;
};

//struct state_t{
//  state_t(int id,int tumor_cn,  int percent,  float het_ratio,bool alpha_identifiable, bool het_identifiable, bool tumor_identifiable, float * sd_rr_alpha, float * sd_rr_loh, float  sd_rr_normal);
//  state_t();
//  void update_rr_mean(float offset);
//  int id;
//  float total_alpha;
//  float het_ratio;
//  //float baf_mean_precompute;
//  bool alpha_ident;
//  bool het_ident;
//  bool tumor_ident;
//  //int alpha_block;
//  int normal_bac; // B allele count
//  int normal_cn;
//  int tumor_bac;
//  int tumor_cn;
//  float mu_rr;
//  float sd_rr;
//  float alpha_prior;
//  void debug();
//     
//};
