#include<boost/property_tree/ptree.hpp>
#include<boost/property_tree/xml_parser.hpp>
#include<list>

using namespace::std;

class Analyzer;

class global_settings_t{
public:
  bool use_db;
  string db_host;
  string db_user;
  string db_pw;
  string db_db;
  string db_table;
  string pedfile;
  string snpfile;
  string genofile;
  string covariatedatafile;
  string covariateselectionfile;
  string selected_analysis;
  list<boost::property_tree::ptree> ptrees;
  void load(const string & filename, Analyzer * & analyzer);
};

