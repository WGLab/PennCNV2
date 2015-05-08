#include<iostream>
#include<fstream>
#include<cstdlib>
#include<map>
#include<set>

#include<vector>
#ifdef USE_DB
#include<mysql.h>
#endif

using namespace std;

typedef map<string,vector<float> > covar_map_t;
typedef vector<vector<float> > covar_vector_t;
//class global_settings_t;


class plink_data_t{
public:
  bool initialized;
  plink_data_t();
  ~plink_data_t();
  bool snp_major;  
  int veclen;
  int totalpersons;
  int totalsnps;
  int * affection;
  char ** genomatrix;
  map<string,int> rs_map;
  vector<string> rs_list;
  map<string,int> subject_map;
  vector<string> subject_list;
};

class IO{
private:
  #ifdef USE_DB
  MYSQL * mysql;
  #endif
  void init();
  void fetchgeno(char * genomatrix,int veclen,float * geno,int totalcols);
  bool plink_format;
  //global_settings_t * settings;
  void loadgarymatrix();
  void loadplinkmatrix();
  //covar_map_t covar_map;

public:
  const static int AFF_CASE = 1;
  const static int AFF_CONTROL = -1;
  const static int AFF_UNDEF = -99;
  const static float MISSING_GENO=9.;
  //const char * prefix;
  virtual ~IO();
  IO();
  //IO(global_settings_t * settings, bool plink_format);
  //IO(int totalpersons,int veclen);
  //IO(const char * prefix); //PLINK style
  void connectToDB(const char * host,const char * user,const char * pw,const char * db);
  bool is_plink_format();

  int snpindex(string rs,plink_data_t & data);
  void readPlinkFiles(const char * bim,const char * fam, const char * bed, plink_data_t & data);
  bool readAllCovariates(const char * data, covar_vector_t & covar_vector);
  bool readCovariates(const char * data,const char * selected, covar_map_t & covar_map);
  void fetchgeno(char * genomatrix,int veclen,int rowindex,float * geno,int totalcols);
  void fetchgeno(char ** genomatrix,int veclen, int rowindex,float * geno,int totalcols);
  bool fetchimputedgeno(const char * db_table,string rsid,float * geno,int len);
  //covar_map_t getcovars();
};
