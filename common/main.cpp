#include<limits.h>
#include<iostream>
#include<fstream>
#ifdef USE_GPU
#include<CL/cl.hpp>
#include"clsafe.h"
#endif
#ifdef USE_MPI
#include<mpi.h>
#endif

//#include<boost/property_tree/exceptions.hpp>
#include"io.hpp"
#include"main.hpp"
#include"analyzer.hpp"
#include"utility.hpp"
#if defined gwas
#include"../gwas/stepwise.hpp"
#include"../gwas/univariate.hpp"
#endif
#if defined pimsa
#include"../pimsa/pathwaysearch.hpp"
#endif
#if defined lasso
#include"../parallel-lasso/dimension2.h"
#include"../parallel-lasso/lasso_mpi2.hpp"
#include"../parallel-lasso/stability.hpp"
#include"../parallel-lasso/power.hpp"
#include"../parallel-lasso/cross_validation.hpp"
#endif
#if defined gpu_impute
#include"../gpu-impute/hmm_impute_dimensions.h"
#include"../gpu-impute/hmm_impute.hpp"
#endif
#if defined cnv
#include"../cnv/hmm_cnv.hpp"
#endif
#if defined ssvs
#include"../ssvs/ssvs_dimensions.h"
#include"../ssvs/ssvs.hpp"
#endif


int main(int argc,char* argv[]){
  if (argc<2){
    cout<<"<analysis>"<<endl;
    exit(0);
  }
  string selected_analysis(argv[1]);
  ostringstream oss;
  oss<<selected_analysis<<".xml";
  string config_file = oss.str();
  ifstream ifs(config_file.data());
  if (!ifs.is_open()){
    cout<<"Configuration file "<<config_file<<" not found.\n";
    exit(0);
  }
  ptree pt;
  read_xml(config_file, pt);
  ifs.close();
  Analyzer * analyzer = NULL;
  //cerr<<"Selected analysis: "<<selected_analysis<<endl;
  try{
    if (!selected_analysis.compare("cross_validation")){
      #if defined lasso
      analyzer = new CrossValidation();
      #endif
    }else if (!selected_analysis.compare("power")){
      #if defined lasso
      analyzer = new Power();
      #endif
    }else if (!selected_analysis.compare("stability")){
      #if defined lasso
      analyzer = new Stability();
      #endif
    }else if (!selected_analysis.compare("univariate")){
      #if defined gwas
      analyzer = new Univariate();
      #endif
    }else if (!selected_analysis.compare("stepwise")){
      #if defined gwas
      analyzer = new Stepwise();
      #endif
    }else if (!selected_analysis.compare("pimsa")){
      #if defined pimsa
      analyzer = new MCMCSampler();
      #endif
    }else if (!selected_analysis.compare("tumor_cnv") || !selected_analysis.compare("seq_cnv") ){
      #if defined cnv
      analyzer = new HMM_CNV(selected_analysis.data());
      #endif
    }else if (!selected_analysis.compare("hmm_impute")){
      #if defined gpu_impute
      analyzer = new HMM_impute();
      #endif
    }else if (!selected_analysis.compare("ssvs")){
      #if defined ssvs
      analyzer = new SSVS();
      #endif
    }
  }catch(const char * & mesg){
    cerr<<"Loader aborted with message: "<<mesg<<endl;
  }
  if (analyzer!=NULL){
    try{
      analyzer->init(pt);
      analyzer->run();
    }catch(const char * & mesg){
      cerr<<"Analyzer aborted with message: "<<mesg<<endl;
    }
    delete(analyzer);
  }else{
    cerr<<"Did not find an analyzer. Exiting\n";
  }
  return 0;
}
