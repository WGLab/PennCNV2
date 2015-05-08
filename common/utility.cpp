#include<iostream>
//#include "gwa.hpp"
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_multifit.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<gsl/gsl_statistics_double.h>
//#include "pathwaysearch.hpp"
//#include "db_manager.hpp"
#include "utility.hpp"


using namespace std;

void debugMatrix(const gsl_matrix *m);


void MathUtils::outputMatrix(const double * mat, uint dim1, uint dim2){
    gsl_matrix_const_view newmat = gsl_matrix_const_view_array(mat,dim1,dim2);
    debugMatrix(&newmat.matrix);
}

double MathUtils::corrIndicator(const double * data1,const double * data2,int len,double * termWeights,double termWeightsTotal){
  double sum=0.;
  for(int i=0;i<len;++i){
    if (data1[i]>0 && data2[i]>0){
      sum+=termWeights[i];
    }
  }
  sum/=(1.*termWeightsTotal);
  return sum;
}

double MathUtils::beta(double a,double b){
  return gsl_ran_beta (rng,a,b);
}

double MathUtils::corr(const double * data1,const double * data2,int len){
    double cor2 = 0.;
    double mu1 = mean(data1,len);
    double mu2 = mean(data2,len);
    for(int i=0;i<len;++i){
      cor2+=1.*(data1[i]-mu1)*(data2[i]-mu2);
    }
    cor2/=1.*(len-1)*SD(data1,len,mu1)*SD(data2,len,mu2);
    //double cor= gsl_stats_correlation(data1,1,data2,1,len);
    if(isnan(cor2)||isinf(cor2)) cor2 = 0.;
    //if(isnan(cor)) cor = -1.;
    //cerr<<"Cor: "<<cor<<", Cor2: "<<cor2<<endl;
    return cor2;
}



void MathUtils::standardize(double * data,uint len){
    double meanbar = mean(data,len);
  //cerr<<"Mean bar: "<<meanbar<<endl;
    double sdbar = SD(data,len,meanbar);
  //cerr<<"SD bar: "<<sdbar<<endl;
    for(uint i=0;i<len;++i){
      data[i]-=meanbar;
      data[i]/=sdbar;
    }
}

MathUtils::MathUtils(unsigned long int seed){
    rng = gsl_rng_alloc (gsl_rng_mt19937);   // instantiate the Mersenne Twister RNG.
    gsl_rng_set(rng,seed); // initialize it with the user desired seed
    cerr<<"RNG type of "<<gsl_rng_name (rng)<<" using seed "<<seed<<endl;
//    printf ("first value = %lu\n", gsl_rng_get (rng));
}


bool MathUtils::chDecomp(const double * inputMatrix, double * decompMat,uint dim){
    for(uint i=0;i<dim*dim;++i) decompMat[i] = inputMatrix[i];
    gsl_matrix_const_view A = gsl_matrix_const_view_array(inputMatrix,dim,dim);
    gsl_matrix_view B = gsl_matrix_view_array(decompMat,dim,dim);
    gsl_matrix_memcpy(&B.matrix,&A.matrix);
//    debugMatrix(&B.matrix);
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    int returncode = gsl_linalg_cholesky_decomp(&B.matrix);
    gsl_set_error_handler (old_handler);
    if (returncode!=0){
        //if (returncode==GSL_EDOM) cerr<<"Not positive definite:\n";
        //  debugMatrix(&A.matrix);
//        cerr<<"Cholesky:\n";
//        debugMatrix(&B.matrix);
//        exit(1);
        return false;

    }
    return true;
    debugMatrix(&B.matrix);
}


double MathUtils::SD(const double * data,uint dim){
    double sd = gsl_stats_sd(data,1,dim);
    if (sd==0){
//        sd = -1;
    }
//    cerr<<"Returning "<<sd<<endl;
    return sd;

}

double MathUtils::SD(const double * data,uint dim,double mean){
    double sd = gsl_stats_sd_m(data,1,dim,mean);
    if (sd==0){
        sd = -1;
    }
//    cerr<<"Returning "<<sd<<endl;
    return sd;

}

double MathUtils::mean(const double * data,uint dim){
    return gsl_stats_mean(data,1,dim);
}

MathUtils::~MathUtils(){
    gsl_rng_free (rng);
}

void MathUtils::ranMVNormal(const double * sigmaMat,
 const double * mu, double * randomSample, uint dimension){

    gsl_permutation * p = gsl_permutation_alloc(dimension);
    gsl_matrix_const_view A = gsl_matrix_const_view_array(sigmaMat,dimension,dimension);
    double decomp[dimension*dimension];
    gsl_matrix_view decompMat = gsl_matrix_view_array(decomp,dimension,dimension);
    gsl_matrix_memcpy(&decompMat.matrix,&A.matrix);
    int signum=0;
    gsl_linalg_LU_decomp(&decompMat.matrix,p,&signum);

//    outputMatrix(decomp,dimension,dimension);
    for(uint i=0;i<dimension;++i){
        randomSample[i] = this->stdNormal(1.0);
    }
    int dim = static_cast<int>(dimension);

    for (int i = dim-1;i >=0; --i) {
        randomSample[i] *= gsl_matrix_get(&decompMat.matrix,i,i);
        for (int j = 0; j < i; ++j) {
            randomSample[i] += gsl_matrix_get(&decompMat.matrix,i,j)* randomSample[j];
        }
    }
    // Add means if necessary
    if (mu != NULL) {
        for (uint i = 0; i < dimension; ++i) {
            randomSample[i] += mu[i];
        }
    }
}



double MathUtils::RandomUniform(){
    return gsl_rng_uniform (this->rng);
}


double MathUtils::stdNormal(double sd){
 double x= gsl_ran_gaussian(this->rng, sqrt(sd));
 return x;
   
}


double MathUtils::cdfNormal(const double mean,const double var,
const double x){
  //cerr<<"mean:"<<mean<<"var:"<<var<<"x:"<<x<<endl;
  double pval =  gsl_cdf_gaussian_P(1.*(mean-x),sqrt(var));
  //double chi = pow(x-mean,2)/var;
  //double pval = (1.-gsl_cdf_chisq_P(chi,1));
  return pval;
  //return gsl_cdf_gaussian_P ((x-mean),sqrt(var));
}

void MathUtils::matrixMultiply(bool transpose1,bool transpose2,double * mat1,double * mat2,double * mat3,int dim1y,int dim1x,int dim2y,int dim2x,int dim3y,int dim3x){
  gsl_matrix_view A = gsl_matrix_view_array(mat1,dim1y,dim1x);
  gsl_matrix_view B = gsl_matrix_view_array(mat2,dim2y,dim2x);
  gsl_matrix_view C = gsl_matrix_view_array(mat3,dim3y,dim3x);
  CBLAS_TRANSPOSE_t t1 = transpose1?CblasTrans:CblasNoTrans;
  CBLAS_TRANSPOSE_t t2 = transpose2?CblasTrans:CblasNoTrans;
  gsl_blas_dgemm(t1,t2,1.0,&A.matrix,&B.matrix,0.,&C.matrix);
}

void MathUtils::invert(const double * decompMat,double * invertedMat,uint dim){
    gsl_matrix_const_view decomp = gsl_matrix_const_view_array(decompMat,dim,dim);
    gsl_matrix_view invmat = gsl_matrix_view_array(invertedMat,dim,dim);
    gsl_matrix_set_identity(&invmat.matrix);
    for(uint i=0;i<dim;++i){
        gsl_vector_view x = gsl_matrix_column (&invmat.matrix, i);
        gsl_linalg_cholesky_svx(&decomp.matrix,&x.vector);
    }
}

bool MathUtils::scoreTest(double * aff, double * designMat, double * betas, double * invInfoMatrix, int samplesize, int rank, double & logL, double & chiSq){
  chiSq = 0.;
  logL = 0.;
  double U[rank];
  double infoMatrix[rank*rank];
  for(int i=0;i<rank;++i){
    U[i]=0;
    for(int j=0;j<rank;++j){
      infoMatrix[i*rank+j]=0;
    }
  }
  //ofstream ofs("test");
  for(int i=0;i<samplesize;++i){
    double betaX = 0.;
  //  ofs<<aff[i];
    for(int j=0;j<rank;++j){
      // compute the inner product
      betaX+=betas[j] * designMat[i*rank+j];
   //   ofs<<"\t"<<designMat[i*rank+j];
    }
    //ofs<<endl;
    double pY = exp(betaX);
    pY/=(1+pY);
    if (aff[i]==1) logL += log(pY); else logL += log(1-pY);
    double dev = aff[i] - pY;
    double info0 = pY * (1-pY);
    for(int j=0;j<rank;++j){
      U[j]+=dev * designMat[i*rank+j];
      for(int k=j;k<rank;++k){
        infoMatrix[j*rank+k]+=info0*designMat[i*rank+j]*designMat[i*rank+k];
        if (k>j) infoMatrix[k*rank+j] = infoMatrix[j*rank+k];
      }
    }
  }
  for(int j=0;j<rank;++j){
    for(int k=0;k<rank;++k){
      //if (k) ofs<<"\t";
      //ofs<<infoMatrix[j*rank+k];
    }
    //ofs<<endl;
  }

  //ofs.close();
  //exit(0);
  double infoMatDecomp[rank*rank];
  if (!chDecomp(infoMatrix,infoMatDecomp,rank)){
    return false;
  }
  invert(infoMatDecomp,invInfoMatrix,rank);
  for(int i=0;i<rank;++i){
    for(int j=0;j<rank;++j){
      betas[i]+=U[j]*invInfoMatrix[i*rank+j];
      chiSq+=1.0*U[i]*invInfoMatrix[i*rank+j]*U[j];
    }
  }
  //cerr<<"Chisq: "<<chiSq<<endl;
  return true;
}


float MathUtils::binom_prob(int n,int k,float p){
  //return 1;
        float num=1;
        for(int i=n;i>(n-k);--i){
                num*=i;
        }
        float den=1;
        for(int i=k;i>1;--i){
                den*=i;
        }
        //cerr<<"num:"<<num<<" den: "<<den<<endl;
        //cerr<<"n: "<<n<<"k: "<<k<<endl;

        float binom = num/den*pow(p,k)*pow(1.-p,n-k);
        //cerr<<"choose: "<<1.*num/den<<" binom: "<<binom<<endl;
        return binom;
}

float MathUtils::gamma_pdf(float x,float shape,float scale){
  if (shape==0) return 0;
  return gsl_ran_gamma_pdf(x,shape,scale);
}

float MathUtils::normal_pdf(float x,float mu,float sigma){
        //float c = 2.506628;
        float p = 1/(sqrt(6.283185*sigma*sigma) * exp(.5 * pow((x-mu)/sigma,2)));
        return p;
}

float MathUtils::normal_cdf (float x, float mu, float sigma){
/* Returns the cumulative probability density function for a normal distribution with mean as mu and standard deviation as sigma
cumulative normal distribution
                  x    2
         1      /   -t  / 2
     ---------- |  e         dt
     sqrt(2 pi) /
                 -inf
*/
   return (1+erf((x-mu)/(sigma*sqrt(2))))/2;
}


bool MathUtils::linReg(double * aff, double * designMat, double * betas, double * invInfoMatrix, int samplesize, int rank, double & logL, double * fitted){
  double xprimex[rank*rank];
  matrixMultiply(true,false,designMat,designMat,xprimex,samplesize,rank,samplesize,rank,rank,rank);
  //cerr<<endl;
  for(int i=0;i<samplesize;++i){
    for(int j=0;j<rank;++j){
    //  if (j) cerr<<" ";
    //  cerr<<designMat[i*rank+j];
    }
   // cerr<<endl;
  }
  double chol[rank*rank];
  if (!chDecomp(xprimex, chol, rank)) return false;
  invert(chol,xprimex,rank); // now xprimex is inverted
  double xxx[rank*samplesize];
  matrixMultiply(false,true,xprimex,designMat,xxx,rank,rank,samplesize,rank,rank,samplesize);
  matrixMultiply(false,false,xxx,aff,betas,rank,samplesize,samplesize,1,rank,1);
  matrixMultiply(false,false,designMat,betas,fitted,samplesize,rank,rank,1,samplesize,1);
  double resid[samplesize];
  double sssr;
  for(int i=0;i<samplesize;++i) resid[i] = aff[i] - fitted[i];
  matrixMultiply(true,false,resid,resid,& sssr,samplesize,1,samplesize,1,1,1);
  double sigma2 = 1.0*sssr/samplesize;
  //sssr/=samplesize;
  for(int i=0;i<rank;++i){
    for(int j=0;j<rank;++j){
      invInfoMatrix[i*rank+j] = sigma2 * xprimex[i*rank+j];
    }
  }
  if (!chDecomp(invInfoMatrix, chol, rank)) return false;
  double invcov[rank * rank];
  invert(chol,invcov,rank); // now xprimex is inverted
  const double ln2pi = 1.837877;
  logL = -0.5*samplesize*ln2pi-0.5*samplesize*log(sigma2)-0.5*samplesize;
  return true;
}

void debugMatrix(const gsl_matrix *m){
    cerr<<"Matrix start\n";
    for(uint i=0;i<m->size1;++i){
        for (uint j=0;j<m->size2;++j){
            if (j) cerr<<" ";
            cerr<<gsl_matrix_get(m,i,j);
        }
        cerr<<endl;
    }
    cerr<<"Matrix end\n";
}



