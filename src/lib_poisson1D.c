/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=-1.0;
    AB[kk+ *kv+1]=2.0;
    AB[kk+ *kv+2]=-1.0;
  }
  AB[0]=0.0;
  if (*kv == 1) {AB[1]=0;}
  
  AB[(*lab)*(*la)-1]=0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}


double relative_forward_error(double* x, double* y, int* la){
  double norm_diff = 0.0;
  double norm_y = 0.0;
  
  // Calculate ||x - y||
  for(int i = 0; i < *la; i++){
      norm_diff += (x[i] - y[i]) * (x[i] - y[i]);
  }
  norm_diff = sqrt(norm_diff);
  
  // Calculate ||y||
  for(int i = 0; i < *la; i++){
      norm_y += y[i] * y[i];
  }
  norm_y = sqrt(norm_y);
  
  // Return ||x - y|| / ||y||
  return norm_diff / norm_y;
}

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
}

int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // Implementation of LU factorization for tridiagonal matrix
  int i, j;
  double m;
  
  // Check input parameters
  if(*kl != 1 || *ku != 1) {
      *info = -1;
      return *info;
  }
  
  // Main loop for LU factorization
  for(i = 0; i < *la-1; i++) {
      // Check for zero pivot
      if(fabs(AB[*lab * i + 1]) < 1e-10) {
          *info = i+1;
          return *info;
      }
      
      // Compute multiplier
      m = AB[*lab * i + 2] / AB[*lab * i + 1];
      AB[*lab * i + 2] = m;
      
      // Update next diagonal element
      AB[*lab * (i+1) + 1] = AB[*lab * (i+1) + 1] - m * AB[*lab * i + 0];
  }
  
  *info = 0;
  return *info;
}
