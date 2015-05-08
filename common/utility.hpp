#include<cmath>
#include<gsl/gsl_randist.h>


class MathUtils{
public:
    MathUtils(unsigned long int seed);
    ~MathUtils();
    double beta(double a,double b);
    bool chDecomp(const double * inputMatrix, double * decompMat,uint dim);
    void invert(const double * decompMat,double * invertedMat,uint dim);
    double SD(const double * data,uint dim);
    double SD(const double * data,uint dim,double mean);
    double corr(const double * data1,const double * data2,int len);
    double corrIndicator(const double * data1,const double * data2,int len,double * termWeights,double termWeightsTotal);
    double mean(const double * data,uint dim);
    void ranMVNormal(const double * sigmaMat, const double * means, double * randomSample, uint dimension);
    void standardize(double * data,uint len);
    void matrixMultiply(bool transpose1,bool transpose2,double * mat1,double * mat2,double * mat3,int dim1y,int dim1x,int dim2y,int dim2x,int dim3y,int dim3x);
    double RandomUniform();
    //double expRV(double dLambda);
    double stdNormal(double sd);
    double cdfNormal(const double mean,const double var,const double x);
    void outputMatrix(const double * mat, uint dim1, uint dim2);
    bool linReg(double * outcome, double * covariates, double * coeff, double * variances, int totalrows, int rank, double & logL, double * fitted);
    bool scoreTest(double * outcome, double * covariates, double * coeff, double * variances, int totalrows, int rank, double & logL, double & chiSq);
    float binom_prob(int n,int k,float p);
    float gamma_function(float k);
    float gamma_pdf(float x,float shape,float scale);
    float normal_pdf(float x,float mu,float sigma);
    float  normal_cdf (float x, float mu, float sigma);
    gsl_rng * rng;
private:
};
