#include<set>
#include<vector>
#include<map>
#include<math.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<gsl/gsl_randist.h>

#ifdef USE_GPU
#include<CL/cl.hpp>
#include"../common/clsafe.h"
#endif
#include"hmm_cnv_dimension.h"
#include"../common/main.hpp"
#include"../common/analyzer.hpp"
#include"../common/io.hpp"
#include"../common/utility.hpp"
#include "hmm_cnv.hpp"
//#include "compute_alpha.hpp"

using namespace std;


void HMM_CNV::init(const ptree & pt){
  verbose_output = false;
  //verbose_output = false;
  LOG_SMALL = log(SMALL);
  LOG_LARGE = log(LARGE);
  this->math = new MathUtils(123);
  // load configurations
  //total_markers = pt.get<int>("total_markers");
  // the following are tuning parameters for sensitivity/specificity
  sd_rr_loh = new float[5];
  sd_baf_percent = new float[100]; 
  sd_rr_alpha = new float[100];
  memset(sd_rr_loh,0,sizeof(float)*5);
  memset(sd_baf_percent,0,sizeof(float)*100);
  memset(sd_rr_alpha,0,sizeof(float)*100);
  if (analysis_type==ANALYSIS_TUMOR){
    this->alpha_percents.push_back(1);
    this->alpha_percents.push_back(33);
    this->alpha_percents.push_back(66);
    this->alpha_percents.push_back(99);
  }else if (analysis_type==ANALYSIS_SEQ){
    this->alpha_percents.push_back(1);
  }
  // 4 states incorporate heterogeneity, normal and LOH state do not.
  int abberant_states = analysis_type==ANALYSIS_TUMOR?4:2;
  this->states = alpha_percents.size() * abberant_states + 2; 
  if(analysis_type==ANALYSIS_TUMOR){
    //
    // TUMOR RELATED TUNING PARAMETERS
    //
    for(vector<int>::iterator it = alpha_percents.begin();
    it!=alpha_percents.end();it++){
      int percent = *it;
      ostringstream oss_tumor;
      oss_tumor<<"lrr_sd_cn.stromal."<<percent;
      sd_rr_alpha[percent]=pt.get<float>(oss_tumor.str());
    }
    for(int i=2;i<=2;++i){
      ostringstream oss;
      oss<<"lrr_sd_cn.loh."<<i;
      sd_rr_loh[i]=pt.get<float>(oss.str());
    }
    sd_rr_normal = pt.get<float>("lrr_sd_normal");
    sd_alpha=pt.get<float>("alpha_sd");
  }else if(analysis_type==ANALYSIS_SEQ){

    //
    // SEQUENCING RELATED TUNING PARAMETERS
    //
    shape = pt.get<float>("seq_info.shape");
    coverage_offset = pt.get<float>("seq_info.mean_coverage_offset");
    dispersion = pt.get<float>("seq_info.dispersion");
  }

  //
  // SHARED TUNING PARAMETERS
  //
  //
  //
  sd_baf_homo_del = pt.get<float>("baf_sd_homo_del");
  normal_cdf = math->normal_cdf(0,.5,sd_baf_homo_del);
  int baf_percents[]={0,25,33,50};
  for(uint i=0;i<sizeof(baf_percents)/sizeof(float);++i){
    ostringstream oss;
    oss<<"baf_sd_percent."<<baf_percents[i];
    sd_baf_percent[baf_percents[i]]=pt.get<float>(oss.str());
  }

  output_dir=pt.get<string>("output_dir");
  cerr<<"Will save results in "<<output_dir<<endl;
  pfb_filename=pt.get<string>("pfb_file");
  cerr<<"Reading population B allele frequencies from "<<pfb_filename<<endl;
  // OK now parse the PFB file if it exists
  ifstream ifs_pfb(pfb_filename.data());
  if (!ifs_pfb.is_open()){
    cerr<<"Warning: cannot locate a population B allele frequency file.  Will assume .5 for missing positions\n";
  }else{
    string line;
    getline(ifs_pfb,line); // skip header
    while(getline(ifs_pfb,line)){
      istringstream iss(line);
      string name,chr;
      int pos;
      float pfb;
      iss>>name>>chr>>pos>>pfb;
      pfb_map[pos] =  pfb;
    }
    ifs_pfb.close();
  }
  input_filelist = pt.get<string>("input_filelist");
  cerr<<"Loading signal files defined in "<<input_filelist<<"\n";
  if (!load_input_filelist(input_filelist)){
    cerr<<"Could not load the set of input files\n";
    exit(1);
  }
}

void HMM_CNV::free_sample_matrices(){
  delete[] this->backtrace_vec;
  delete[] this->pos_vec;
  //this->alpha_block_vec = new int[total_markers];
  delete[] this->rslist;
  delete[] this->gap_vec;
  delete[] this->forward_rescalings;
  delete[] this->backward_rescalings;
  delete[] this->greedy_rescalings;
  delete[] this->intensity_vec;
  delete[] this->segment_length_vec;
  delete[] this->baf_vec;
  delete[] this->emission_cache_mat;
  delete[] this->forward_mat;
  delete[] this->greedy_mat;
  delete[] this->bestpath_mat;
  delete[] this->backward_mat;
  delete[] this->prob_mat;
  delete[] this->prev_prob_mat;
}

void HMM_CNV::allocate_global_matrices(){
  cerr<<"Total States: "<<states<<endl;
  state_het_ratio = new float[states];
  state_total_alpha = new float[states];
  state_tumorcn_vec = new int[states];
  state_het_ident = new bool[states];
  state_alpha_ident = new bool[states];
  state_tumor_ident =  new bool[states];
  state_alpha_prior = new float[states];
  state_mu_rr = new float[states];
  state_sd_rr = new float[states];
  //sample_likelihood_log = new float[samples];
  trans_matrix_size = states * states;
  this->trans_mat = new float[trans_matrix_size];
}

void HMM_CNV::allocate_sample_matrices(){
  cerr<<"Scanning "<<sample_info_vec[current_sample].filename<<" for number of markers...\n";
  ifstream ifs_test(sample_info_vec[current_sample].filename.data());
  if (!ifs_test.is_open()){
    cerr<<"Could not open "<<sample_info_vec[current_sample].filename<<endl;
    throw "File not found.";
  }
  total_markers = 0;
  string testline;
  // skip header
  getline(ifs_test,testline);
  // count lines after header
  while(getline(ifs_test,testline)) ++total_markers;
  ifs_test.close();
  cerr<<total_markers<<" found.\n";
  this->backtrace_vec = new int[total_markers];
  this->pos_vec = new int[total_markers];
  this->rslist = new string[total_markers];
  this->gap_vec = new int[total_markers];
  this->forward_rescalings = new int[total_markers];
  this->backward_rescalings = new int[total_markers];
  this->greedy_rescalings = new int[total_markers];
  this->intensity_vec = new float[total_markers];
  this->segment_length_vec = new int[total_markers];
  this->baf_vec = new float[total_markers];
  hmm_matrix_size = states * total_markers;
  this->emission_cache_mat = new float[hmm_matrix_size];
  this->forward_mat = new float[hmm_matrix_size];
  this->greedy_mat = new float[hmm_matrix_size];
  this->bestpath_mat = new int[hmm_matrix_size];
  this->backward_mat = new float[hmm_matrix_size];
  this->prob_mat = new float[hmm_matrix_size];
  this->prev_prob_mat = new float[hmm_matrix_size];
};

void HMM_CNV::initialize_state_parameter(int state,int tumor_cn, int percent, float het_ratio, bool alpha_identifiable,bool het_identifiable, bool tumor_identifiable, float * sd_rr_alpha, float * sd_rr_loh, float  sd_rr_normal){
  state_tumorcn_vec[state] = tumor_cn;
  state_het_ident[state] = het_identifiable;
  state_alpha_ident[state] = alpha_identifiable;
  state_tumor_ident[state] =  tumor_identifiable;
  state_het_ratio[state] = het_ratio;
  state_total_alpha[state] = percent/100.;
  float slope = (1.-INTERCEPT)/2.;
  if (intensity_type==INTENSITY_LRR){
    state_mu_rr[state] = (1-state_total_alpha[state])*(INTERCEPT+tumor_cn*slope)+ (state_total_alpha[state]);
    if (tumor_cn==2){
      if (het_ratio>0){
        // the variance of the normal tissue
        state_sd_rr[state] = (sd_rr_normal);
      }else{
        state_sd_rr[state] = (sd_rr_loh[tumor_cn]);
      }
    }else{
      state_sd_rr[state] = sd_rr_alpha[percent];
    }
  }else if (intensity_type==INTENSITY_READS){
    state_mu_rr[state] = tumor_cn/2.;
    state_sd_rr[state] = dispersion;
  }
  cerr<<"for cn "<<tumor_cn<<" with contamination level "<<state_total_alpha[state]<<", mu_rr: "<<state_mu_rr[state]<<", and sd_rr: "<<state_sd_rr[state]<<endl;
  float alpha_prior = intensity_type==INTENSITY_LRR?math->normal_pdf(state_total_alpha[state],global_alpha,sd_alpha):1;
  if(alpha_prior<SMALL) alpha_prior = SMALL;
  state_alpha_prior[state] =  alpha_prior; 
}

//void HMM_CNV::initialize_state_parameters(){
  // because each state has means that are a function a sample
  // specific stromal contamination, these must be initialized
  // on a sample by sample basis.
//  init_state_objects(current_sample);
  // set up normal state
  //for(int state=0;state<states;++state){
//    state_het_ratio[state] = state_object_vector[state].het_ratio;
//    state_total_alpha[state] = state_object_vector[state].total_alpha;
//    state_mu_rr[state] = state_object_vector[state].mu_rr;
//    //cerr<<"state: "<<state<<" has mu_rr "<<state_mu_rr[state]<<endl;
//    state_tumorcn_vec[state] = state_object_vector[state].tumor_cn;
//    state_het_ident[state] = state_object_vector[state].het_ident;
//    state_alpha_ident[state] = state_object_vector[state].alpha_ident;
//    state_sd_rr[state] = state_object_vector[state].sd_rr;
//  }
//}

void HMM_CNV::initialize_global_matrices(){
}

HMM_CNV::HMM_CNV(const char * analysis_type){
  cerr<<"Selected analysis type of "<<analysis_type<<endl;
  if (strcmp(analysis_type,"tumor_cnv")==0){
    this->analysis_type = ANALYSIS_TUMOR;
    this->intensity_type = INTENSITY_LRR;
  }else if (strcmp(analysis_type,"seq_cnv")==0){
    this->analysis_type = ANALYSIS_SEQ;
    this->intensity_type = INTENSITY_READS;
  }else{
    throw "Unknown analysis type\n";
  }
}

bool HMM_CNV::load_input_filelist(const string & filelistpath){
  ifstream ifslist(filelistpath.data());
  if (!ifslist.is_open()){
    cerr<<"File list is not found.\n";
    return false;
  }
  string line;
  // skip header
  getline(ifslist,line);
  float epsilon=.01;
  int sample_id = 0;
  while(getline(ifslist,line)){
    sample_info_t sample_info;
    istringstream iss(line);
    if (intensity_type==INTENSITY_LRR){
      iss>>sample_info.filename>>sample_info.annot>>sample_info.global_alpha;
      if (sample_info.global_alpha==0) sample_info.global_alpha+=epsilon;
      else if (sample_info.global_alpha==1) sample_info.global_alpha-=epsilon;
    }else{
      iss>>sample_info.filename>>sample_info.annot;
    }
    ++sample_id;
    sample_info_vec.push_back(sample_info);
  }
  ifslist.close();
  samples = sample_info_vec.size();
  cerr<<"Total samples: "<<samples<<endl;
  if (samples==0) {
    cerr<<"Must have at least one sample\n";
    return false;
  }
  return true;
}

void HMM_CNV::load_input_file(){
  //map<int,int> existing_position_map;
  ifstream ifsobs(sample_info_vec[current_sample].filename.data());
  int lastpos = 0;
  string line;
  getline(ifsobs,line); // skip header
  string last_chr="";
  int last_pos = -1;
  //cerr<<"Reading in "<<total_markers<<" lines.\n";
  for(int j=0;j<total_markers;++j){
    getline(ifsobs,line);
    istringstream iss(line);
    float intensity;
    switch (intensity_type){
      case INTENSITY_LRR:
        iss>>rslist[j]>>chr>>pos_vec[j]>>baf_vec[j]>>intensity;
        if (last_chr.compare("")!=0 && last_chr.compare(chr)!=0){
          cerr<<"Found chr string of "<<chr<<" and prev chr of "<<last_chr<<endl;
          throw "Error in input: all rows must be from the same chromosome\n";
        }
        last_chr = chr;
        if (last_pos!=-1 && pos_vec[j]<=last_pos){
          cerr<<"Found SNP of "<<rslist[j]<<" and position of "<<pos_vec[j]<<" and previous position of "<<last_pos<<endl;
          throw "Error in input: all rows must contain positions in increasing order\n";
        }
        last_pos = pos_vec[j];
        intensity_vec[j] = pow(2,intensity);
        break;
      case INTENSITY_READS:
        iss>>rslist[j]>>chr>>pos_vec[j]>>baf_vec[j]>>intensity_vec[j]>>segment_length_vec[j];
        break;
    }
    gap_vec[j] = pos_vec[j]-lastpos;
    lastpos = pos_vec[j];
    // for the first sample, record all observed positions
    //if (i==0) {
      //cerr<<"Mapping position "<<pos_vec[j]<<" as index "<<j<<endl;
      //existing_position_map[pos_vec[j]] = j;
    //}
  }
  ifsobs.close();
  cerr<<"Loaded "<<sample_info_vec[current_sample].filename.data()<<endl;
  if(intensity_type ==INTENSITY_READS){
    coverage = get_mean_depth();
    cerr<<"Empirical read depth is "<<coverage<<endl;
    coverage+=coverage_offset;
    cerr<<"New user revised mean depth is "<<coverage<<endl;
  }
}

void HMM_CNV::run(){
  allocate_global_matrices();
  for(current_sample=0;current_sample<samples;++current_sample){
    allocate_sample_matrices();
    load_input_file();
    initialize_state_parameters();
    cerr<<"\n\nBeginning HMM_CNV on sample "<<sample_info_vec[current_sample].annot<<" with stromal contamination of "<<global_alpha<<endl;
    if (intensity_type==INTENSITY_LRR) adjust_aneuploidy();
    compute_emission_cache_mat();
    // before baum welch training, initialize trans mat
    float last_loglike = likelihood_log = -1.*LARGE;
    iteration = 0;
    //for(iteration=0;iteration<5;++iteration){
    bool rollback = false;
    do{
      last_loglike = likelihood_log  ;
      for(int i=0;i<total_markers;++i){
        for(int j=0;j<states;++j){
          prev_prob_mat[i*states+j] = prob_mat[i*states+j];
          //cerr<<" "<<trans_mat[i*states+j];
        }
        //cerr<<endl;
      }
      if (iteration==0) {
        this->forwardbackward();
      }else{
        this->baumwelch();
        this->forwardbackward();
      } 
      cerr<<"Last LL: "<<last_loglike<<" and current: "<<likelihood_log<<endl;
      ++iteration;
      if (likelihood_log<last_loglike){
        rollback = true;
      }
    //}
    //}while(iteration<2);
    }while(!rollback);
    cerr<<"Rolling back to previous probability matrix\n";
    for(int i=0;i<total_markers;++i){
      for(int j=0;j<states;++j){
        prob_mat[i*states+j] = prev_prob_mat[i*states+j] ;
        //cerr<<" "<<trans_mat[i*states+j];
      }
      //cerr<<endl;
    }
    viterbi();
    print_summary();
    free_sample_matrices();
  }
}


HMM_CNV::~HMM_CNV(){
  cerr<<"Destructor\n";
  delete[] sd_rr_alpha;
  delete[] sd_rr_loh;
  delete[] sd_baf_percent;
}

void HMM_CNV::adjust_aneuploidy(){
    float mu_rr_offset = get_rr_offset();
    cerr<<"Adjusted RR mean for normal state is now: "<<(1. + mu_rr_offset) <<endl;
    for(int state=0;state<states;++state){
      state_mu_rr[state]+=mu_rr_offset;
    }
}

float HMM_CNV::get_binom_sum(int n, float pop_b, float baf){
  float binom_sum = 0;
  float sd = .1;
  float mean=.5;
  float mean_increment = 0;
  if (n>1){
  	mean = 0;
  	mean_increment =1./(n-1);
  }
  for(int k=0;k<n;++k){
  	binom_sum+=math->binom_prob(n,k,pop_b) * math->normal_pdf(baf,mean,sd);
  	mean+=mean_increment;
  }
  return binom_sum;
}



float ordered_binom(int n,int k,float p){
  float prob = 1;
  for(int i=0;i<k;++i){
    prob*=p;
  }
  for(int i=0;i<n-k;++i){
    prob*=(1-p);
  }
  return prob;
}

// this needs to be called several times for each state

float HMM_CNV::get_baf_mean(int bac, int cn,float alpha){
  if((cn>=2 && bac==0) || bac==cn){
    // if tumor homo , assume normal is always homo
    return bac/cn;
  }
  int normal_bac = 1;
  int normal_cn = 2;
  float mu = ((1.-global_alpha-alpha)*bac + (global_alpha+alpha)*normal_bac)
      / ((1.-global_alpha-alpha)*cn+(global_alpha+alpha)*normal_cn);
  //cerr<<"alpha "<<alpha<<" tumor bac,cn "<<bac<<","<<cn<<" bafmu: "<<mu<<endl;
  return mu;
}


float HMM_CNV::get_baf_sd(int state, int bac){
  //float & alpha = state_alpha[state];
  float & total_alpha = state_total_alpha[state];
  // the variance of the normal tissue
  int tumorcn = state_tumorcn_vec[state];
  //cerr<<"\nFetched tumor CN "<<tumorcn<<" with het "<<state_het_ratio[state]<<" and BAC "<<bac<<endl;
  float frac = 1.*bac/tumorcn;
  //cerr<<"frac "<<frac<<endl;
  if (frac>.5) frac = 1.-frac;
  int percent = (int)(frac * 100. + .5);
  float sdnorm = sd_baf_percent[50];
  float sdtumor = sd_baf_percent[percent];
  float var = pow((total_alpha)*sdnorm,2)
   + pow((1-total_alpha)*sdtumor,2);
  //cerr<<"debug: "<<pow((1-global_alpha-alpha)*sdtumor,2)<<endl;
  //cerr<<"alpha "<<global_alpha<<","<<alpha<<" SDs "<<sdnorm<<","<<sdtumor<<"; variance: "<<var<<endl;
  return sqrt(var);
}


//state_t::state_t(int id,int tumor_cn, int percent, float het_ratio, bool alpha_identifiable,bool het_identifiable, bool tumor_identifiable, float * sd_rr_alpha, float * sd_rr_loh, float  sd_rr_normal){
//  this->id = id;
//  this->tumor_cn = tumor_cn;
//  this->het_ident = het_identifiable;
//  this->alpha_ident = alpha_identifiable;
//  this->tumor_ident = tumor_identifiable;
//  this->het_ratio = het_ratio;
//  this->total_alpha = percent/100.;
//  float slope = (1.-INTERCEPT)/2.;
//  this->mu_rr = (1-total_alpha)*(INTERCEPT+tumor_cn*slope)+ (total_alpha);
//  //cerr<<"for cn "<<tumor_cn<<" mu_rr: "<<this->mu_rr<<endl;
//  if (tumor_cn==2){
//    if (het_ratio>0){
//      // the variance of the normal tissue
//      sd_rr = (sd_rr_normal);
//    }else{
//      sd_rr = (sd_rr_loh[tumor_cn]);
//    }
//  }else{
//    sd_rr = sd_rr_alpha[percent];
//  }
//}

//void state_t::update_rr_mean(float offset){
//  this->mu_rr += offset;
//}
//
//void state_t::debug(){
//  cerr<<id<<"\t"<<tumor_cn<<"\t"<<het_ratio<<"\t"<<total_alpha<<"\t"<<alpha_ident<<het_ident<<tumor_ident<<"\t"<<mu_rr<<"\t"<<sd_rr<<"\t"<<alpha_prior;
//  cerr<<endl;
//}

void HMM_CNV::initialize_state_parameters(){
  normal_state = -1;
  homo_deletion_state = -1;
  //state_object_vector.clear();
  global_alpha = intensity_type==INTENSITY_LRR?sample_info_vec[current_sample].global_alpha:0;
  int states = 0 ;
  int percent = global_alpha*100;
  // normal state
  initialize_state_parameter(states,2,percent,HET_RATIO,false,true,true, sd_rr_alpha,sd_rr_loh,sd_rr_normal);
  //state_object_vector.push_back(state_t(states,2,percent,HET_RATIO,false,true,true, sd_rr_alpha,sd_rr_loh,sd_rr_normal));
  if (normal_state==-1) normal_state = states;
  ++states;
  // copy neutral LOH
  initialize_state_parameter(states,2,percent,0,false,true, false,sd_rr_alpha,sd_rr_loh,sd_rr_normal);
  ++states;
  for(vector<int>::iterator it = alpha_percents.begin();
  it!=alpha_percents.end();it++){
    int percent = *it;
    if (analysis_type==ANALYSIS_TUMOR){
      // homozygous deletion
      initialize_state_parameter(states,0,percent,HET_RATIO,true,false, true,sd_rr_alpha,sd_rr_loh,sd_rr_normal);
      ++states;
    }
    // hemizygous deletion
    initialize_state_parameter(states,1,percent,HET_RATIO,true,false, true,sd_rr_alpha,sd_rr_loh,sd_rr_normal);
    ++states;
    // amplifications
    initialize_state_parameter(states,3,percent,HET_RATIO,true,false, true,sd_rr_alpha,sd_rr_loh,sd_rr_normal);
    ++states;
    if (analysis_type==ANALYSIS_TUMOR){
      initialize_state_parameter(states,4,percent,HET_RATIO,true,false, true,sd_rr_alpha,sd_rr_loh,sd_rr_normal);
      ++states;
    }
  }
  assert(states == this->states);
  cerr<<"Total states matched: "<<states<<endl;
  for(int i=0;i<states;++i){
    for(int j=0;j<states;++j){
      trans_mat[i*states+j] = 1./(states);
    }
  }
  cerr<<"Transition matrix initialized to uninformative matrix\n";
}

void HMM_CNV::normalize(int len, float * raw_probs, float * normalized){
  float normalizer = 0;
  for (int j=0;j<len;++j){
    normalizer+=raw_probs[j];
  }
  for (int j=0;j<len;++j){
    normalized[j] = raw_probs[j]/normalizer;
    if (normalized[j]<SMALL) normalized[j] = SMALL;
  }
}

void HMM_CNV::normalize_across_alphas(float * baf_emission_vec){
  for(int cn = 0;cn<TUMOR_CN_ELEMENTS;++cn){
    if (cn!=2){
      float normalizer = 0;
      for(int state = 0;state<states;++state){
        if (state_tumorcn_vec[state]==cn){
           normalizer+=baf_emission_vec[state];
        }
      }
      for(int state = 0;state<states;++state){
        if (state_tumorcn_vec[state]==cn){
           baf_emission_vec[state]/=normalizer;
           //cerr<<"new baf for state "<<state<<": "<<baf_emission_vec[state]<<endl;
        }
      }
    }else{
      // this is a uniform probability across all alpha percents
      baf_emission_vec[normal_state] = 1./alpha_percents.size();
    }
  }
}

void HMM_CNV::compute_emission_cache_mat(){
  bool debug = false;
  cerr<<"Computing and storing emission probabilities\n";
  //ofstream ofs_emission("debugging/emission_cache");
  ostringstream oss;
  oss<<output_dir<<"/emission.debug."<<current_sample;
  ofstream ofs_debug;
  if (debug) ofs_debug.open(oss.str().data());
  for (int i=0;i<total_markers;++i){
    float pfb = pfb_map.find(pos_vec[i])==pfb_map.end()?.5:
    pfb_map[pos_vec[i]];
    if (verbose_output && i%10000==0) cerr<<" "<<i;
    float & intensity = intensity_vec[i];
    float & baf = baf_vec[i];
    //float baf_noise = .01;
    float intensity_noise= .01;
    float probs[states];
    float baf_probs[states];
    float intensity_probs[states];
    for (int state=0;state<states;++state){
      if (debug) ofs_debug<<"Marker "<<i<<" at position "<<pos_vec[i]<<",";
      // get emission of the tumor heterogeneity
      float & alpha = state_total_alpha[state];
      //float & alpha_prob = state_alpha_prior[state];
      // figure out the emission probability of the LRR
      int tumorcn = state_tumorcn_vec[state];
      float & intensity_prob = intensity_probs[state];
      if (intensity_type==INTENSITY_LRR){
        intensity_prob = math->normal_pdf(intensity,state_mu_rr[state],state_sd_rr[state]);
      }else if (intensity_type==INTENSITY_READS){
        // set this to one to ignore length scaling 
        // for now.  
        float scaling_coeff = 1;
        //float scaling_coeff = segment_length_vec[i]; 
        float scale = state_sd_rr[state]/scaling_coeff;
        // the expected value of a gamma distrib is shape * scale
        float mean = shape * scale;
        // how far is the observed value from expected value for the state
        float deviance = intensity - (coverage * state_mu_rr[state]) ;
        // add offset of mean as appropriate for the gamma dist
        float expected_intensity = mean+deviance;
        intensity_prob = math->gamma_pdf(expected_intensity,shape,scale);
        //cerr<<"at marker "<<i<<" pos "<<pos_vec[i]<<" state "<<state<<" intensity: "<<intensity<<" mean "<<mean<<" deviance "<<deviance<<" exp_intensity "<<expected_intensity<<" shape "<<shape<<" scale "<<scale<<endl;
      }
      intensity_prob = intensity_noise+(1.-intensity_noise)*intensity_prob;
      //cerr<<"at marker "<<i<<" state "<<state<<" intensity: "<<intensity<<" prob "<<intensity_prob<<endl;
      //cerr<<"Final SD and prob for state "<<state<<": "<<state_sd_lrr[state]*sd_lrr_inflation<<","<<lrr_prob<<endl;
      // figure out the emission prob of the BAF    
      float & baf_prob = baf_probs[state];
      baf_prob = 0;
      float pointmass=.5;
      if (tumorcn==0){
        if (baf==0||baf==1){
          baf_prob+=normal_cdf;
        }else{
          baf_prob+=math->normal_pdf(baf,.5,sd_baf_homo_del);
        }
      }else{  // these other cases take into account stromal contamination
        bool status_loh = (tumorcn==1 || (tumorcn==2 && state_het_ratio[state]==0));
        int effective_tumorcn = status_loh?1:tumorcn;
        // the first two cases only happen when tumor and normal are homozygote
        if (baf==0){
          baf_prob+=math->binom_prob(effective_tumorcn,0,pfb)*pointmass;
        }else if (baf==1){
          baf_prob+=math->binom_prob(effective_tumorcn,effective_tumorcn,pfb)
          * pointmass;
        }else{
          // this loop assumes the normal genotype is a heterozygote.
          for (int bac=0;bac<=effective_tumorcn;++bac){
            baf_prob+=math->binom_prob(effective_tumorcn,bac,pfb)*math->normal_pdf(baf,get_baf_mean(bac,effective_tumorcn,alpha),get_baf_sd(state,bac));
          }
          if (status_loh){
            // assume homozygote normal too
            int bac = (baf<.5)?0:1;
            baf_prob+=math->binom_prob(effective_tumorcn,bac,pfb)*math->normal_pdf(baf,bac,sd_baf_percent[0]);
          }
        }
      }
      //cerr<<"observation: "<<i<<endl;
      //baf_prob = baf_noise+(1.-baf_noise)*baf_prob;
      //probs[state] = baf_prob;
      //probs[state] = rr_prob * alpha_prob;
      //probs[state] = rr_prob * baf_prob * alpha_prob;
      if (baf_prob<SMALL) baf_prob = SMALL;
      if (intensity_prob<SMALL) intensity_prob = SMALL;
      if (debug){
        ofs_debug<<"ID:"<<state<<","<<"BAF:"<<baf<<",intensity:"<<intensity<<","<<"P_Intensity:"<<intensity_prob<<","<<"P_BAF:"<<baf_prob<<endl;
      }
    }
    normalize_across_alphas(baf_probs);
    for(int state=0;state<states;++state){
      probs[state] = intensity_probs[state]*baf_probs[state]*state_alpha_prior[state];
      //cerr<<"ID:"<<state<<","<<"BAF:"<<baf<<",intensity:"<<intensity<<","<<"P_Intensity:"<<intensity_probs[state]<<","<<"P_BAF:"<<baf_probs[state]<<"alpha: "<<state_alpha_prior[state]<<endl;
      //emission_cache_mat[i*states+state] = rr_probs[state]*baf_probs[state]*state_alpha_prior[state];
    }
    normalize(states,probs,emission_cache_mat+i*states);
  }
  if (verbose_output) cerr<<endl;
  //ofs_emission.close();
  if (debug) ofs_debug.close();
  for(int j=0;j<states;++j){
    forward_mat[j] = emission_cache_mat[j];
    greedy_mat[j] = emission_cache_mat[j];
    backward_mat[(total_markers-1)*states+j] = 1;
  }
  cerr<<"Done initializing emission_cache_mat\n";
}


double HMM_CNV::addMultipleLogs(double * logs, int len){
  //cerr<<"Called addMultipleLogs\n";
  double avg=0;
  for(int i=0;i<len;++i){
    avg-=logs[i];
  }
  avg/=len;
  //avg-=10;
  //cerr<<"Average scaling: "<<avg<<endl;
  double val = 0;
  for(int i=0;i<len;++i){
    val+=exp(logs[i]+avg);
  }
  double a = log(val)-avg;
  //cerr<<"Results "<<a<<endl;
  return a;
}


void HMM_CNV::print_summary(){
  ostringstream oss1;
  ostringstream oss2;
  ostringstream oss3;
  ostringstream oss4;
  ostringstream oss5;
  ostringstream oss6;
  oss1<<output_dir<<"/summary."<<sample_info_vec[current_sample].annot;
  oss2<<output_dir<<"/posterior_alpha."<<sample_info_vec[current_sample].annot;
  oss3<<output_dir<<"/posterior_het."<<sample_info_vec[current_sample].annot;
  oss4<<output_dir<<"/posterior_tumor_cn."<<sample_info_vec[current_sample].annot;
  oss5<<output_dir<<"/transition."<<sample_info_vec[current_sample].annot;
  oss6<<output_dir<<"/abberation."<<sample_info_vec[current_sample].annot;


  cerr<<"Marginalizing for sample "<<sample_info_vec[current_sample].annot<<" across "<<total_markers<<" markers"<<endl;
  ofstream ofs_expected(oss1.str().data());
  ofstream ofs_alpha(oss2.str().data());
  ofstream ofs_het(oss3.str().data());
  ofstream ofs_tumor_cn(oss4.str().data());
  ofstream ofs_trans(oss5.str().data());
  ofstream ofs_abb(oss6.str().data());
  for(int i=0;i<states;++i){
    ofs_trans<<i;
    for(int j=0;j<states;++j){
      ofs_trans<<"\t"<<j<<":"<<trans_mat[i*states+j];
    }
    ofs_trans<<endl;
  }
  ofs_trans.close();
  float alpha_probs[100];
  float het_probs[100];
  float tumor_cn_probs[TUMOR_CN_ELEMENTS];
  if (intensity_type==INTENSITY_LRR){
    ofs_expected<<"chr\tpos\tBAF\tRR\ttumor_cn(expected)\ttumor_cn(best)\talpha\theterozyg\n";
    ofs_abb<<"chr\tpos\tlength\ttumor_cn(best)\n";
  }else if (intensity_type==INTENSITY_READS){
    ofs_expected<<"chr\tpos\tBAF\tRR\tcn(expected)\tcn(best)\theterozyg\n";
    ofs_abb<<"chr\tpos\tlength\tcn(best)\n";
  }
  int last_cn = state_tumorcn_vec[backtrace_vec[0]];
  int start_seg = pos_vec[0];
  int end_seg = pos_vec[1];
  
  for(int i = 0;i<total_markers;++i){
    ofs_alpha<<i;
    ofs_het<<i;
    ofs_tumor_cn<<i;
    memset(alpha_probs,0,sizeof(alpha_probs));
    memset(het_probs,0,sizeof(het_probs));
    memset(tumor_cn_probs,0,sizeof(tumor_cn_probs));
    for(int state=0;state<states;++state){
      //state_t s = state_object_vector[state];
      float prob = prob_mat[i*states+state];
      //cerr<<"state "<<state<<" id "<<s.alpha_ident<<" alpha: "<<s.total_alpha<<endl;
      if(state_alpha_ident[state]){
        int alpha_int = (int)(state_total_alpha[state]*100);
        alpha_probs[alpha_int]+=prob;
      }
      if(state_het_ident[state]){
        int het_int = (int)(state_het_ratio[state]*100);
        het_probs[het_int]+=prob;
      }
      if(state_tumor_ident[state]){
        // if normal state, scale posterior by number of fraction states to be comparable to abberation posteriors
        tumor_cn_probs[state_tumorcn_vec[state]]+=(state==normal_state)?prob*alpha_percents.size():prob;
      }
    }
    float alpha_normalizer = 0;
    float het_normalizer = 0;
    for(int j=0;j<100;++j){
      if (alpha_probs[j]>0) alpha_normalizer += alpha_probs[j];
      if (het_probs[j]>0) het_normalizer += het_probs[j];
    }
    float exp_alpha = 0;
    float exp_het = 0;
    for(int j=0;j<100;++j){
      if (alpha_probs[j]>0){
         alpha_probs[j]/=alpha_normalizer;
         ofs_alpha<<"\t"<<(j/100.)<<":"<<alpha_probs[j];
         exp_alpha+=alpha_probs[j] * (j/100.);
      }
      if (het_probs[j]>0){
         het_probs[j]/=het_normalizer;
         ofs_het<<"\t"<<(j/100.)<<":"<<het_probs[j];
         exp_het+=het_probs[j] * (j/100.);
      }
    }
    float normalizer = 0;
    for(int j=0;j<TUMOR_CN_ELEMENTS;++j){
      if (tumor_cn_probs[j]>0) normalizer += tumor_cn_probs[j];
    }
    float exp_tumor_cn=0;
    for(int j=0;j<TUMOR_CN_ELEMENTS;++j){
      if (tumor_cn_probs[j]>0) {
        tumor_cn_probs[j]/=normalizer;
        ofs_tumor_cn<<"\t"<<j<<":"<<tumor_cn_probs[j];
        exp_tumor_cn+=tumor_cn_probs[j] * j;
      }
    }
    ofs_alpha<<endl;
    ofs_het<<endl;
    ofs_tumor_cn<<endl;
    int best_cn = state_tumorcn_vec[backtrace_vec[i]];
    if (intensity_type==INTENSITY_LRR){
      ofs_expected<<chr<<"\t"<<pos_vec[i]<<"\t"<<baf_vec[i]<<"\t"<<(intensity_vec[i])<<"\t"<<exp_tumor_cn<<"\t"<<best_cn<<"\t"<<exp_alpha<<"\t"<<exp_het<<endl;
    }else if (intensity_type==INTENSITY_READS){
      ofs_expected<<chr<<"\t"<<pos_vec[i]<<"\t"<<baf_vec[i]<<"\t"<<(intensity_vec[i])<<"\t"<<exp_tumor_cn<<"\t"<<best_cn<<"\t"<<exp_het<<endl;
    }
    if (best_cn==last_cn){
      end_seg = i<(total_markers-1)?pos_vec[i+1]:pos_vec[i];
    }else{
      int length = end_seg-start_seg;
      ofs_abb<<chr<<"\t"<<start_seg<<"\t"<<length<<"\t"<<last_cn<<endl;
      last_cn = best_cn;
      start_seg =  pos_vec[i];
      end_seg = i<(total_markers-1)?pos_vec[i+1]:pos_vec[i];
    }
  }
  int length = end_seg-start_seg;
  ofs_abb<<chr<<"\t"<<start_seg<<"\t"<<length<<"\t"<<last_cn<<endl;
  ofs_alpha.close();
  ofs_het.close();
  ofs_tumor_cn.close();
  ofs_expected.close();
  ofs_abb.close();
}

void HMM_CNV::forwardbackward(){
  cerr<<"Running forward backward...\n";
  for(int i=0;i<total_markers;++i){
    forward_rescalings[i] = 0;
    backward_rescalings[i] = 0;
  }
  //ofstream ofs_forward("debugging/forward");
  // RUN FORWARD ALGORITHM
  for(int obs=1;obs<total_markers;++obs){
    forward_rescalings[obs] = forward_rescalings[obs-1];
  //cerr<<"Launching obs"<<total_markers<<"\n";
    if (verbose_output && obs % 10000 == 0) cerr<<" "<<obs;
    float min = 1;
    for(int currentstate=0;currentstate<states;++currentstate){
      float sum = 0;
      for(int prevstate=0;prevstate<states;++prevstate){
         sum+=forward_mat[(obs-1)*states+prevstate]  * trans_mat[prevstate*states+currentstate] * emission_cache_mat[obs*states+currentstate];
      }
      if (sum<min) min = sum;
      forward_mat[obs*states+currentstate] = sum;
    }
    if (min<SMALL){
      for(int currentstate=0;currentstate<states;++currentstate){
        forward_mat[obs*states+currentstate]*=LARGE;
      }
      //if (verbose_output) cerr<<"Rescaled forward probabilities at marker "<<obs<<endl;
      ++forward_rescalings[obs];
    }
  }
  // initialize the probabilities at the last observation as 1 since no data beyond
  // RUN BACKWARD ALGORITHM
  for(int obs=total_markers-2;obs>=0;--obs){
    backward_rescalings[obs] = backward_rescalings[obs+1];
    if (verbose_output && obs % 10000 == 0) cerr<<" "<<obs;
    float min = 1;
    for(int currentstate=0;currentstate<states;++currentstate){
      float sum = 0;
      for(int futurestate=0;futurestate<states;++futurestate){
        sum+= backward_mat[(obs+1)*states+futurestate] * trans_mat[currentstate*states+futurestate]  * emission_cache_mat[(obs+1)*states+futurestate];
      }
      if (sum<min) min = sum;
      backward_mat[obs*states+currentstate] = sum;
    }
    if (min<SMALL){
      for(int currentstate=0;currentstate<states;++currentstate){
        backward_mat[obs*states+currentstate]*=LARGE;
      }
      //if (verbose_output) cerr<<"Rescaled backward probabilities at marker "<<obs<<endl;
      ++backward_rescalings[obs];
    }
  }
  if (verbose_output) cerr<<endl;
  if (verbose_output) cerr<<"Forward and backward rescalings: "<<forward_rescalings[total_markers-1]<<","<<backward_rescalings[0]<<endl;
  // COMPUTE LOG LIKELIHOODS
  double sum_log[states];
  double sum_log2[states];
  for(int state=0;state<states;++state){
    sum_log[state]=log(forward_mat[state])+log(backward_mat[state])+backward_rescalings[0]*LOG_SMALL;
    sum_log2[state]=log(forward_mat[(total_markers-1)*states+state])+log(backward_mat[(total_markers-1)*states+state])+forward_rescalings[total_markers-1]*LOG_SMALL;
  }
  float likelihood_log = addMultipleLogs(sum_log,states);
  float likelihood_log2 = addMultipleLogs(sum_log2,states);
  cerr<<"Iteration: "<<iteration<<": CPU likelihood_logs: "<<likelihood_log<<","<<likelihood_log2<<endl;
  if (abs(likelihood_log-likelihood_log2)>1){
    cerr<<"Mismatch in likelihoods. See debugging dir\n";
    ofstream ofs_emission("debugging/emission.mat");
    ofstream ofs("debugging/posterior.mat");
    ofstream ofs_f("debugging/forward.mat");
    ofstream ofs_b("debugging/backward.mat");
    for(int obs=0;obs<total_markers;++obs){
      ofs<<obs;
      ofs_f<<obs;
      ofs_b<<obs;
      ofs_emission<<obs<<"\t"<<"Intensity:"<<(intensity_vec[obs])<<"\t"<<"BAF:"<<baf_vec[obs];
      for(int i=0;i<states;++i){
        ofs<<"\t"<<prob_mat[obs*states+i];
        ofs_f<<"\t"<<forward_mat[obs*states+i];
        ofs_b<<"\t"<<backward_mat[obs*states+i];
        ofs_emission<<"\t"<<i<<":"<<emission_cache_mat[obs*states+i];
      }
      ofs<<endl;
      ofs_f<<endl;
      ofs_b<<endl;
      ofs_emission<<endl; 
    }
    ofs.close();
    ofs_f.close();
    ofs_b.close();
    ofs_emission.close();
    exit(1);
  }
  this->likelihood_log = likelihood_log;
  for(int obs=0;obs<total_markers;++obs){
    for(int state=0;state<states;++state){
      prob_mat[obs*states+state] = exp(log(forward_mat[obs*states+state])+forward_rescalings[obs]*LOG_SMALL+log(backward_mat[obs*states+state])+backward_rescalings[obs]*LOG_SMALL-likelihood_log); 
    }
  }
  cerr<<"Completed computation of posterior probability matrix\n";
  bool debug2 = false;
  if (debug2){
    ofstream ofs_prob("debugging/posterior.mat2");
    ofs_prob<<"CPU Probability matrix: \n";
    for(int obs=0;obs<total_markers;++obs){
      //float rowsum=0;
      ofs_prob<<"obs:"<<obs;
      for(int state=0;state<states;++state){
        ofs_prob<<" "<<prob_mat[obs*states+state];
        //rowsum+=exp(prob_mat_log[obs*states+state]);
      }
    ofs_prob<<endl;
    //cerr<<" "<<rowsum<<endl;
    }
    ofs_prob.close();
    exit(0);
  }
}

void HMM_CNV::baumwelch(){
  cerr<<"Running Baum Welch"<<endl;
  float denom_f_arr[states];
  for(int i=0;i<states;++i){
    denom_f_arr[i] = 0;
    for(int obs=0;obs<total_markers-1;++obs){
      denom_f_arr[i]+=prob_mat[obs*states+i];
    }
  }
  bool debug1 = false;
  if (debug1){
    cerr<<"CPU denom:\n";
    for(int i=0;i<states;++i){
      cerr<<" "<<denom_f_arr[i];
    }
    cerr<<endl;
    //exit(0);
  }
  // new transition kernel
  float numer_f[states*states];
  for(int i=0;i<states*states;++i) numer_f[i]=0.;
  //float likelihood_log=likelihood_log;
  bool debug1a = false;
  ofstream ofs_debug;
  if (debug1a) ofs_debug.open("debugging/baum.txt");
  int max_rescalings = forward_rescalings[total_markers-1] + backward_rescalings[0];
  bool defined[max_rescalings];
  float normalizer_arr[max_rescalings];
  for(int i=0;i<max_rescalings;++i) defined[i] = false;
  for(int obs=0;obs<total_markers-1;++obs){
    if (verbose_output && obs % 1000 == 0) cerr<<" "<<obs;
    float normalizer = 0;
    int index = forward_rescalings[obs]+backward_rescalings[obs];
    if (defined[index]){
      normalizer = normalizer_arr[index];
    }else{
      normalizer = normalizer_arr[index] = exp(forward_rescalings[obs]*LOG_SMALL + backward_rescalings[obs]*LOG_SMALL - likelihood_log);
      defined[index] = true;
    }
    for(int i=0;i<states;++i){
      for(int j=0;j<states;++j){
        float val1 = 0,val2 = 0;
        val2=forward_mat[obs*states+i]*trans_mat[i*states+j] * emission_cache_mat[(obs+1)*states+j] * backward_mat[(obs+1) * states + j] * normalizer;
     
        if(debug1a) ofs_debug<<obs<<"\t"<<i<<"\t"<<j<<"\t"<<val1<<"\t"<<val2<<endl;
        numer_f[i*states+j]+=val2;
      }
    }
  }
  if (debug1a) ofs_debug.close();
  if (verbose_output) cerr<<endl;
  for(int i=0;i<states;++i){
    for(int j=0;j<states;++j){
      //cerr<<"numer "<<i<<","<<j<<": "<<numer_f[i*states+j]<<endl;
      trans_mat[i*states+j]=(numer_f[i*states+j]/denom_f_arr[i]);
      if (trans_mat[i*states+j]<SMALL) trans_mat[i*states+j] = SMALL;
      //trans_mat[i*states+j]= HMM_SMALL + (1.-HMM_SMALL)*trans_mat[i*states+j];
    }
  }
  
  bool debug2 = false;
  if (debug2){
    cerr<<"CPU transition\n";
    for(int i=0;i<states;++i){
      float rowsum = 0;
      for(int j=0;j<states;++j){
        rowsum+=trans_mat[i*states+j];
        cerr<<" "<<trans_mat[i*states+j];
      }
      cerr<<" rowsum: "<<rowsum<<endl;
    }
  }
}

float HMM_CNV::get_mean_depth(){
  float depth = 0;
  for(int i=0;i<total_markers;++i){
    depth += intensity_vec[i];
  }
  depth/=1.*total_markers;
  return depth;
}

float HMM_CNV::get_rr_offset(){
  //return 0;
  float rr_offset = 0;
  float denom = 0;
  float epsilon = .01;
  for(int i=0;i<total_markers;++i){
    float & rr = intensity_vec[i];
    float & baf = baf_vec[i];
    //float rr = pow(2,lrr);
    if (abs(baf-.5)<epsilon){
      float weight = math->normal_pdf(baf,.5,.001); 
      rr_offset+=weight * (rr-1.);
      denom+=weight; 
    }
  }
  rr_offset = denom>0?rr_offset/denom:0;
  return rr_offset;
}

void HMM_CNV::viterbi(){
  cerr<<"Running Viterbi decoding...\n";
  //float sum,forwardsum,backwardsum;
  //int back_trace[total_markers][states];
  //int trace[total_markers];
  float max = SMALL;
  int obs = 0;
  for(int j=0;j<states;++j){
     //cerr<<" "<<greedy_mat[j];
     if (greedy_mat[j] > max){
       backtrace_vec[obs] = j;
       max = greedy_mat[j];
     }
  }
  for(int i=0;i<total_markers;++i) greedy_rescalings[i] = 0;
  // RUN VITERBI
  //ofstream ofs_g("debugging/greedy.mat");
  for(int obs=1;obs<total_markers;++obs){
    //ofs_g<<obs;
    greedy_rescalings[obs] = greedy_rescalings[obs-1];
    if (verbose_output && obs % 10000 == 0) cerr<<" "<<obs;
    float min = 1;
    for(int currentstate=0;currentstate<states;++currentstate){
      float max = 0;
      int maxind = -1;
      for(int prevstate=0;prevstate<states;++prevstate){
        //cerr<<"prevstate "<<prevstate<<" currstate "<<currentstate<<" Got here\n";
        //double curr_val = greedy_mat[(obs-1)*states+prevstate] ;
        float curr_val = greedy_mat[(obs-1)*states+prevstate] * trans_mat[prevstate*states+currentstate];
        //cerr<<"prevstate "<<prevstate<<" currstate "<<currentstate<<" Got here\n";
        if (curr_val>=max){
          maxind = prevstate;
          max = curr_val;
        }
   
      }
      if (maxind<0) {
        cerr <<"Never found a max beyond "<<max<<endl;
        exit(1);
      }
      if (max<min) min = max;
      greedy_mat[obs*states+currentstate] = max * emission_cache_mat[obs*states+currentstate];
      //ofs_g<<"\t"<<greedy_mat[obs*states+currentstate]<<","<<emission_cache_mat[obs*states+currentstate]<<","<<maxind<<endl;
      //cerr<<obs<<" "<<currentstate<<" Max state "<<maxind<<" max prob: "<<max<<" greedy "<<greedy_mat[obs*states+currentstate]  <<" emission: "<<log(emission_cache_mat[obs*states+currentstate]) <<endl;
      bestpath_mat[obs*states+currentstate] = maxind;
    }
    //ofs_g<<endl;
    if (min<SMALL){
      for(int currentstate=0;currentstate<states;++currentstate){
        greedy_mat[obs*states+currentstate]*=LARGE;
      }
      //if (verbose_output) cerr<<"Rescaled forward probabilities at marker "<<obs<<endl;
      ++greedy_rescalings[obs];
    }
  }
  //ofs_g.close();
  backtrace_vec[total_markers-1] = 0;
  float currval = SMALL;
  for(int state=0;state<states;++state){
    if (greedy_mat[(total_markers-1)*states+state] > currval){
      backtrace_vec[total_markers-1] = state;
      currval = greedy_mat[(total_markers-1)*states+state];
    }
  }
  for(int obs = total_markers-2;obs>=0;--obs){
    backtrace_vec[obs] = bestpath_mat[(obs+1)*states+backtrace_vec[obs+1]];
  }
}

