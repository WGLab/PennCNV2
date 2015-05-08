#include<string.h>

#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<boost/property_tree/exceptions.hpp>

#include"analyzer.hpp"
#include"io.hpp"

using namespace std;

const int MAX=100;

Analyzer::Analyzer(){ }

//void Analyzer::run(){}

Analyzer::~Analyzer(){
  cerr<<"Analyzer is deleted.\n";
}

bool Analyzer::chDecomp(const double * inputMatrix, double * decompMat,uint dim){
    for(uint i=0;i<dim*dim;++i) decompMat[i] = inputMatrix[i];
    gsl_matrix_const_view A = gsl_matrix_const_view_array(inputMatrix,dim,dim);
    gsl_matrix_view B = gsl_matrix_view_array(decompMat,dim,dim);
    gsl_matrix_memcpy(&B.matrix,&A.matrix);
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    int returncode = gsl_linalg_cholesky_decomp(&B.matrix);
    gsl_set_error_handler (old_handler);
    if (returncode!=0){
        return false;
    }
    return true;
}

void Analyzer::invert(const double * decompMat,double * invertedMat,uint dim){
    gsl_matrix_const_view decomp = gsl_matrix_const_view_array(decompMat,dim,dim);
    gsl_matrix_view invmat = gsl_matrix_view_array(invertedMat,dim,dim);
    gsl_matrix_set_identity(&invmat.matrix);
    for(uint i=0;i<dim;++i){
        gsl_vector_view x = gsl_matrix_column (&invmat.matrix, i);
        gsl_linalg_cholesky_svx(&decomp.matrix,&x.vector);
    }
}

void Analyzer::outputmat(double * mat,int dim1,int dim2){
  for(int i=0;i<dim1;++i){
    for(int j=0;j<dim2;++j){
      if (j) cerr<<" ";
      cerr<<mat[i*dim2+j];
    }
    cerr<<endl;
  }
  exit(0);
}

void Analyzer::matrixMultiply(bool transpose1,bool transpose2,double * mat1,double * mat2,double * mat3,int dim1y,int dim1x,int dim2y,int dim2x,int dim3y,int dim3x){
  gsl_matrix_view A = gsl_matrix_view_array(mat1,dim1y,dim1x);
  gsl_matrix_view B = gsl_matrix_view_array(mat2,dim2y,dim2x);
  gsl_matrix_view C = gsl_matrix_view_array(mat3,dim3y,dim3x);
  CBLAS_TRANSPOSE_t t1 = transpose1?CblasTrans:CblasNoTrans;
  CBLAS_TRANSPOSE_t t2 = transpose2?CblasTrans:CblasNoTrans;
  gsl_blas_dgemm(t1,t2,1.0,&A.matrix,&B.matrix,0.,&C.matrix);
}

bool Analyzer::leastsquares(double* designMat, double * aff, double * betas, double * cov,double * resid,int rank,int samplesize){
//void matrixMultiply(bool transpose1,bool transpose2,double * mat1,double * mat2,double * mat3,int dim1y,int dim1x,int dim2y,int dim2x,int dim3y,int dim3x)
  double xprimex[rank*rank];
  matrixMultiply(true,false,designMat,designMat,xprimex,samplesize,rank,samplesize,rank,rank,rank);

  double chol[rank*rank];
  if (!chDecomp(xprimex, chol, rank)) return false;
  invert(chol,xprimex,rank); // now xprimex is inverted
  double xxx[rank*samplesize];
  matrixMultiply(false,true,xprimex,designMat,xxx,rank,rank,samplesize,rank,rank,samplesize);
  matrixMultiply(false,false,xxx,aff,betas,rank,samplesize,samplesize,1,rank,1);
  double fitted[samplesize];
  matrixMultiply(false,false,designMat,betas,fitted,samplesize,rank,rank,1,samplesize,1);
  for(int i=0;i<samplesize;++i) resid[i] = aff[i] - fitted[i];
  double sssr;
  matrixMultiply(true,false,resid,resid,& sssr,samplesize,1,samplesize,1,1,1);
  sssr/=(samplesize-rank);
  for(int i=0;i<rank;++i) cov[i*rank+i] = sssr * xprimex[i*rank+i];
  //cerr<<sqrt(cov[rank*rank-1])<<endl; exit(0);
  return true;

}



void Analyzer::scoreTest(double & chiSq,double & logL,int * diseaseStatus,
float * effSizeMat,double * betas,double * invInfoMatrix,uint iObsSampleSize,
suint stride,suint params, int & count){
    chiSq = logL = 0.;
    double U[params];
    double infoMatrix[params*params];
    memset(U,0,sizeof(U));
    memset(infoMatrix,0,sizeof(infoMatrix));
    //ofstream ofs("debug");
    for(uint i=0;i<iObsSampleSize;++i){
        double betaX = 0.;
        uint col = i*stride;
        //ofs<<diseaseStatus[i];
        for(uint j=0;j<params;++j){
           //ofs<<" "<<effSizeMat[col+j];
            betaX+=betas[j] * effSizeMat[col+j];
          //cerr<<"betaX"<<betaX<<endl;
        }
        //ofs<<endl;
        double pY = exp(betaX);
        pY/=(1+pY);
        if (diseaseStatus[i]==1) logL += log(pY); else logL += log(1-pY);
        double dev = diseaseStatus[i] - pY;
        double info0 = pY * (1-pY);
        
        for(uint j=0;j<params;++j){
            U[j]+=dev * effSizeMat[col+j];
            for(uint k=j;k<params;++k){
                infoMatrix[j*params+k]+=info0*effSizeMat[col+j]*effSizeMat[col+k];
                //cerr<<" "<<infoMatrix[j*params+k];
                if (k>j) infoMatrix[k*params+j] = infoMatrix[j*params+k];
            }
            //cerr<<endl;
        }
    }
    //ofs.close();
    //exit(1);
    double infoMatDecomp[params*params];
    if (!chDecomp(infoMatrix,infoMatDecomp,params)){
        cerr<<"Cholesky failed\n";
        ofstream ofs1("dataset");
        ofs1<<"aff";
        for (uint j=0;j<params;++j){
          ofs1<<" p"<<j;
        }
        ofs1<<endl;
        for(uint i=0;i<iObsSampleSize;++i){
          ofs1<<diseaseStatus[i];
          for (uint j=0;j<params;++j){
            ofs1<<" "<<effSizeMat[i*stride+j];
          }
          ofs1<<endl;
        }
        ofs1.close();
        ofstream ofs2("infomat");
        for (uint i=0;i<params;++i){
          for (uint j=0;j<params;++j){
            if (j) ofs2<<" ";
            ofs2<<infoMatrix[i*params+j];
          }
          ofs2<<endl;
        }
        ofs2.close();
        //exit(0);
        count = MAX;
        return;
    }
    invert(infoMatDecomp,invInfoMatrix,params);
    for(uint i=0;i<params;++i){
        uint col = i*params;
        for(uint j=0;j<params;++j){
            betas[i]+=U[j]*invInfoMatrix[col+j];
            chiSq+=1.0*U[i]*invInfoMatrix[col+j]*U[j];
        }
    }
}

bool Analyzer::fitModelCCD(double & logL, int * aff, float * X,double * beta, double * beta_var, int n,int stride, int rank){
  float lambda =0;
  int iter=0;
  double tolerance = 0.;
  double delta_beta[rank];
  double delta_score[n];
  double score[n];
  float X_t[n*stride];
  float Y[n];
  for (int i=0;i<n;++i){
    Y[i] = aff[i]==1?1:-1;
  }
  for (int i=0;i<n;++i){
    for (int j=0;j<rank;++j){
      X_t[j*n+i] = X[i*stride+j];
    }
  }
  for (int j=0;j<rank;++j){
    beta[j] = 0.;
    delta_beta[j] = 0.;
  }
  for (int i=0;i<n;++i){
    score[i] = 0.;
    delta_score[i] = 0.;
  }
  do{
    for (int j=0;j<rank;++j){
      float num = 0;
      float den = 0;
      for (int i=0;i<n;++i){
        float escore = exp(score[i]);
        num += Y[i]*X_t[j*n+i]/(1+escore);
        den += X_t[j*n+i]*X_t[j*n+i]*escore/((1+escore)*(1+escore));
      }
      if(beta[j]==0.){
        delta_beta[j] = (num-n*lambda)/den;
        if (delta_beta[j]<0){
          delta_beta[j] = (num+n*lambda)/den;
          if (delta_beta[j]>0){
            delta_beta[j] = 0;
          }
        }
      }else{
        delta_beta[j] = (num-beta[j]/fabs(beta[j])*n*lambda)/den;
        if ((beta[j]>0 && beta[j]+delta_beta[j]<0) ||
            (beta[j]<0 && beta[j]+delta_beta[j]>0 )){
          delta_beta[j] = -1*beta[j];
        }
      }
      if (delta_beta[j]!=0){
        beta[j]+=delta_beta[j];
        for(int i=0;i<n;++i){
          delta_score[i] = delta_beta[j]*X_t[j*n+i]*Y[i];
          score[i]+=delta_score[i];
        }
      }
    } // end loop through variables
    tolerance = 0;
    double totalcorr = 0;
    for(int i=0;i<n;++i){
      tolerance+=fabs(delta_score[i]);
      totalcorr+=fabs(score[i]);
    }
    for(int j=0;j<rank;++j){
      cout<<" "<<beta[j];
    }
    cout<<endl;
    //cerr<<"Iter: "<<iter<<", Tolerance: "<<tolerance<<", Corr: "<<totalcorr<<endl;
    tolerance/=(1.+totalcorr);
    ++iter;
  }while(tolerance>.000001);
  cerr<<"Iterations: "<<iter;
  return true;
}

bool Analyzer::fitModel(double & L1,int * phenovec_filtered,float * designmat,double * betas, double * var, int samplesize,int stride, int rank){
  int count=0;
  double chiSq;
  //cerr<<"Rank: "<<rank<<" stride: "<<stride<<endl;
  do{
    scoreTest(chiSq,L1,phenovec_filtered,designmat,betas,var,samplesize,stride,rank,count);
    //cerr<<"Iteration: "<<count<<" ChiSq: "<<chiSq<<" logL "<<L1<<endl;
  }while(chiSq>.0001 && count++<MAX);
  if (count>=MAX) {
    return false;
  }else{
    return true;
  }
}

