#include<iostream>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include"io.hpp"
#include"main.hpp"

using namespace std;

//int masks[4] = {3,12,48,192};
//int totalpersons;
//int totalsnps;

//int veclen;
//char ** genomatrix;
const int PLINK_WT_HOMO=0;
const int PLINK_MUT_HOMO=3;
const int PLINK_HET=2;
const int PLINK_MISSING=1;

plink_data_t::plink_data_t(){
  initialized = false;
  affection = NULL;
  genomatrix = NULL;
}

plink_data_t::~plink_data_t(){
  if (affection!=NULL) delete[]affection;
  int totalrows = snp_major?totalsnps:totalpersons;
  if (genomatrix!=NULL){
    for(int i=0;i<totalrows;++i){
      delete[]genomatrix[i]; 
    }
    delete[]genomatrix;
  }
  cerr<<"PLINK data object deleted\n";
}
IO::IO(){}
IO::~IO(){cerr<<"IO object deleted\n";}

#ifdef USE_DB
void IO::connectToDB(const char * host,const char * user,const char * pw,const char * db){
  mysql=NULL;
  if (mysql_library_init(0, NULL, NULL)) {
    throw "Could not initialize MySQL library";
  }
  mysql = mysql_init(mysql);
  if (mysql==NULL){
    throw "MySQL cannot init";
  }
  cerr<<"Using host "<<host<<endl;
  mysql = mysql_real_connect(mysql,host,user,pw,db,3306,NULL,0);
  if (mysql==NULL){
    throw "MySQL cannot connect";
  }
}
#endif


int IO::snpindex(string rsid, plink_data_t & data){
  map<string,int>::iterator it = data.rs_map.find(rsid);
  if(it==data.rs_map.end()) return -1;
  return it->second;
}

bool IO::readAllCovariates(const char * data, covar_vector_t & covar_vector){
  // NOW LOAD THE RAW DATA
  ifstream ifs_covdata(data);
  if (ifs_covdata.is_open()){
    string line,token;
    getline(ifs_covdata,line);
    istringstream iss(line);
    while(iss>>token){
      vector<float> subvec;
      covar_vector.push_back(subvec);
    }
    int totalcols = covar_vector.size();
    while(getline(ifs_covdata,line)){
      istringstream iss(line);
      float cell;
      for(int i=0;i<totalcols;++i){
        iss>>cell;
        covar_vector[i].push_back(cell);
      }
    }
    ifs_covdata.close();
  }
  cerr<<"Covariates loaded\n";
  return true;
  // let's run a test!
  //covar_map_t::iterator it = covar_map.find("e1");
//  vector<float> vec = it->second;
//  for(uint i=0;i<vec.size();++i){
//    cerr<<" "<<vec[i];
//  }
//  
//   exit(1);
}

bool IO::readCovariates(const char * data,const char * selected, covar_map_t & covar_map){
  if (data==NULL || selected==NULL){
    cerr<<"No covariate filenames specified.  Not loading\n";
    return true;
  }
  // NOW LOAD SELECTED PHENOTYPES
  set<string> sel_cov_set;
  ifstream ifs_seldata(selected);
  string line;
  if (ifs_seldata.is_open()){
    while(getline(ifs_seldata,line)){
      sel_cov_set.insert(line);
    }
    ifs_seldata.close();
  }else{
    cerr<<"Found a covariate data file but not selection file\n";
    return false;
  }
  cerr<<"Selecting "<<sel_cov_set.size()<<" covariates\n";
  // NOW LOAD THE RAW DATA
  ifstream ifs_covdata(data);
  if (ifs_covdata.is_open()){
    string header,colname;
    getline(ifs_covdata,line);
    istringstream iss(line);
    vector<string> namevec;
    while(iss>>colname){
      namevec.push_back(colname);
    }
    int totalcols = namevec.size();
    bool include[totalcols];
    for(int i=0;i<totalcols;++i){
      set<string>::iterator it = sel_cov_set.find(namevec[i]);
      include[i] = it!=sel_cov_set.end()?true:false;
      if (include[i]) {
        vector<float> values;
        covar_map[namevec[i]] = values;
      }
    }
    while(getline(ifs_covdata,line)){
      istringstream iss(line);
      string cell;
      for(int i=0;i<totalcols;++i){
        iss>>cell;
        if (include[i]){
          covar_map_t::iterator it = covar_map.find(namevec[i]);
          (it->second).push_back(atof(cell.data()));
        }
      }
    }
    ifs_covdata.close();
  }
  cerr<<"Covariates loaded\n";
  return true;
  // let's run a test!
  //covar_map_t::iterator it = covar_map.find("e1");
//  vector<float> vec = it->second;
//  for(uint i=0;i<vec.size();++i){
//    cerr<<" "<<vec[i];
//  }
//  
//   exit(1);
}

void IO::readPlinkFiles(const char * bim,const char * fam, const char * bed, plink_data_t & data ){
  if (!data.initialized && (bim==NULL||fam==NULL||bed==NULL)){
    cerr<<"Data was not initialized and PLINK file names missing, not loading\n";
    throw "Reading PLINK files failed";
  }
  string line;
  if (fam!=NULL){
    data.totalpersons = 0;
    ifstream ifsperson(fam);
    while(getline(ifsperson,line)){
      istringstream iss(line);
      string subject_id;
      iss>>subject_id;
      data.subject_map[subject_id] = data.totalpersons;
      data.subject_list.push_back(subject_id);
      ++data.totalpersons;
    } 
    cerr<<"Found "<<data.totalpersons<<" persons.\n";
    ifsperson.close();
    // read case control status
    data.affection = new int[data.totalpersons];
    ifsperson.open(fam);
    string val1;
    int j = 0;
    int aff;
    while(getline(ifsperson,line)){
      istringstream iss(line);
      for(uint i=0;i<5;++i) iss>>val1;
      iss>>aff;
      if (aff==1){
        data.affection[j] = AFF_CONTROL;
      }else if (aff==2){
        data.affection[j] = AFF_CASE;
      }else{
        data.affection[j] = AFF_UNDEF;
      }
      ++j;
    } 
    cerr<<"Loaded "<<(j)<<" affection status rows\n";
    ifsperson.close();
  }
  if (bim!=NULL){
    data.totalsnps = 0;
    ifstream ifssnp(bim);
    while(getline(ifssnp,line)){
      istringstream iss(line);
      int chr,position;
      string rsid;
      float morgan;
      char other,ref;
      iss>>chr>>rsid>>morgan>>position>>other>>ref;
      data.rs_map[rsid] = data.totalsnps;
      data.rs_list.push_back(rsid);
      ++data.totalsnps;
    }
    cerr<<"Found "<<data.totalsnps<<" snps.\n";
    ifssnp.close();
  }
  if (bed!=NULL){
    // NOW READ GENOTYPES
    ifstream ifs(bed);
    if (!ifs.is_open()){
       cerr<<"Cannot find "<<bed<<endl;
       exit(0);
    }
    // read header for SNP major mode
    char header[3];
    ifs.read(header,3);
    data.snp_major = (int)header[2];
    cerr<<"SNP major?: "<<data.snp_major<<endl;
    int rows = data.snp_major?data.totalsnps:data.totalpersons;
    int cols = data.snp_major?data.totalpersons:data.totalsnps;
  
    bool remainder = (cols % 4 !=0)?true:false;
    data.veclen = cols/4+remainder; // for remainder
    char vector_read[data.veclen];
    //cerr<<"Allocating "<<data.veclen<<"\n";
    data.genomatrix = new char * [rows];
    for (int row=0;row<rows;++row){
      data.genomatrix[row] = new char[data.veclen];
      ifs.read(vector_read,data.veclen);
      //cerr<<" "<<row<<endl;
      for (int col=0;col<data.veclen;++col){
        data.genomatrix[row][col] = vector_read[col];
      }
    }
    ifs.close();
    cerr<<"Loaded genotype bit matrix of "<<rows<<" by "<<cols<<endl;
  }
  data.initialized = true;
}

void IO::init(){
}

//IO::IO(int totalpersons,int veclen){
//  init();
//  this->totalpersons = totalpersons;
//  this->veclen = veclen;
//}


//void IO::loadgarymatrix(){
//  ostringstream bedstr;
//  //bedstr<<this->prefix<<".gary";
//  ifstream ifs(bedstr.str().data());
//  int totalinlen = totalpersons+1; // for EOF
//  char inputvec[totalinlen];
//  
//  for (int snp=0;snp<totalsnps;++snp){
//    //bool remainder = (totalpersons % 4 !=0)?true:false;
//    //int totaloutlen = totalpersons/4+remainder; // for remainder
//    char  * outputvec = genomatrix[snp];
//    ifs.read(inputvec,totalinlen);
//    int outputcount=0;
//    int totalshifts = 0;
//    for (int person=0;person<totalpersons;++person){
//      if (totalshifts==0) outputvec[outputcount] = 0;
//      int plinkgeno;
//      switch(inputvec[person]){
//        case '1':
//          plinkgeno = PLINK_WT_HOMO;
//          break;
//        case '2':
//          plinkgeno = PLINK_HET;
//          break;
//        case '3':
//          plinkgeno = PLINK_MUT_HOMO;
//          break;
//        case '0':
//          plinkgeno = PLINK_MISSING;
//          break;
//      }
//      //cerr<<"plinkgeno: "<<plinkgeno;
//      //cerr<<"person "<<person<<" snp "<<snp<<" geno "<<inputvec[snp]<<endl;
//      //cerr<<"t_shifts "<<totalshifts<<" outputcount "<<outputcount<<" plinkgeno: "<<plinkgeno;
//      //cerr<<" b_shift"<<plinkgeno;
//      for(int shift=0;shift<totalshifts;++shift) plinkgeno = plinkgeno<<2;
//      //cerr<<" a_shift"<<plinkgeno;
//      outputvec[outputcount] = outputvec[outputcount] | plinkgeno;
//      //cerr<<" after add: "<<(int)outputvec[outputcount];
//      if (totalshifts==3){
//        totalshifts = 0;
//        ++outputcount;
//      }else{
//        ++totalshifts;
//      }
//    }
//    //cout.write(outputvec,totaloutlen);
//    //ofs.write(outputvec,totaloutlen);
//    //cerr<<"Person "<<person<<" completed.\n";
//  }
//  ifs.close();
//  //ofs.close();
//}


#ifdef USE_DB
bool IO::fetchimputedgeno(const char * db_table,string rsid,float * geno,int len){
  // attempt a look up in the database
  ostringstream select1;
  select1<<"select vector from "<<db_table
  <<" where rs='"<<rsid<<"'";
  //cerr<<select1.str().data()<<endl;
  mysql_query(mysql,select1.str().data());
  MYSQL_RES * res = mysql_store_result(mysql);
  if (res){
    MYSQL_ROW row = mysql_fetch_row(res);
    int count=0;
    if (row){
      string str=string(row[count++]);
      for(int i=0;i<len;++i){
        geno[i] = ((int)str[i]-65)/10.;
      }
    }else{
      cerr<<"Can not locate rsid"<<rsid<<endl;
      return false;
    }
  }
  mysql_free_result(res);
  return true;
}
#endif

//void IO::fetchgeno(int snpindex,float * geno){
//  fetchgeno(this->genomatrix,snpindex,geno);
//}

void IO::fetchgeno(char * genomatrix,int veclen, int rowindex,float * geno, int len){
  unsigned int i = rowindex*veclen;
  char * vector_read = genomatrix+i;
  fetchgeno(vector_read,veclen,geno,len);
}

void IO::fetchgeno(char ** genomatrix,int veclen, int rowindex,float * geno,int len){
  char * vector_read = genomatrix[rowindex];
  fetchgeno(vector_read,veclen,geno,len);
}


// used internally, shared by all the above function signatures.

void IO::fetchgeno(char * vector_read,int veclen,float * geno,int totalcols){
  int masks[4];
  masks[0] = 3;
  masks[1] = 12;
  masks[2] = 48;
  masks[3] = 192;
  int colindex = 0;
  //cerr<<"Vec len: "<<veclen<<" total persons: "<<totalcols<<endl;
  for (int byteindex =0;byteindex<veclen;++byteindex){
    for(int pair=0;pair<4;++pair) {
      int val = masks[pair] & vector_read[byteindex];
      for(int shift=0;shift<pair;++shift) val = val>>2;
      if (colindex<totalcols)  {
        switch (val){
        case PLINK_WT_HOMO:
          geno[colindex] = 0;
          break;
        case PLINK_HET:
          geno[colindex] = 1;
          break;
        case PLINK_MUT_HOMO:
          geno[colindex] = 2;
          break;
        case PLINK_MISSING:
          geno[colindex] = MISSING_GENO;
          break;
        default:
          cerr<<"Not recognized: "<<geno[colindex]<<endl;
          exit(0);
          break;
     
        }
        //cerr<<" "<<val<<","<<geno[colindex];
      } 
      ++colindex;
    }
  }
}

//covar_map_t IO::getcovars(){
//  return covar_map;
//}



int main2(int argc,char * argv[]){
  if (argc<2){
    cerr<<"Arguments: <fileprefix>\n";
    exit(1);
  }
  //int totalpersons = atoi(argv[1]);
  //int totalsnps = atoi(argv[2]);
  //IO io(NULL,true);
  //io.readInputFiles();
  //io.loadplinkmatrix();
  //float test[io.totalpersons];
  //io.fetchgeno(1,test);
  //for(int i=0;i<io.totalpersons;++i){
   //cout<<test[i];
  //}
  //cout<<endl;
  return 0;
}
