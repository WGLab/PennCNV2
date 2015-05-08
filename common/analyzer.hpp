class IO;

typedef unsigned int uint;
typedef short unsigned int suint;

// provides several convenience functions for all analyzers
using boost::property_tree::ptree;

class Analyzer{
public:
  Analyzer();
  virtual void init(const ptree & pt)=0;
  virtual void run()=0;
  virtual ~Analyzer();
protected:
  bool fitModelCCD(double & logLikelihood,int * affection_status,float * designmat,double * betas, double * beta_covar, int samplesize,int stride, int rank);
  bool fitModel(double & logLikelihood,int * affection_status,float * designmat,double * betas, double * beta_covar, int samplesize,int stride, int rank);
  const char * config_file;
  void matrixMultiply(bool transpose1,bool transpose2,double * mat1,double * mat2,double * mat3,int dim1y,int dim1x,int dim2y,int dim2x,int dim3y,int dim3x);
  void outputmat(double * mat,int dim1,int dim2);
  bool leastsquares(double * designMat, double * aff, double * betas, double * cov,double * resid,int rank,int samplesize);
  bool chDecomp(const double * inputmat,double * decompmat,uint dim);
  void invert(const double * d,double * i,uint dim);
  void scoreTest(double & chiSq,double & logL,int * diseaseStatus,
float * designmat,double * betas,double * invInfoMatrix,uint iObsSampleSize,
suint stride,suint params, int & count);
};
