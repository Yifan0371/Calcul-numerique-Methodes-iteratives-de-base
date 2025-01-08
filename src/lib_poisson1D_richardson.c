/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

// Calculate eigenvalues for 1D Poisson problem
void eig_poisson1D(double* eigval, int *la) {
  int i;
  double h = 1.0 / ((*la) + 1);
  for (i = 0; i < *la; i++) {
      eigval[i] = 4 * pow(sin((i + 1) * M_PI * h / 2.0), 2);
  }
}

// Calculate maximum eigenvalue
double eigmax_poisson1D(int *la) {
  double h = 1.0 / ((*la) + 1);
  return 4 * pow(sin((*la) * M_PI * h / 2.0), 2);
}

// Calculate minimum eigenvalue
double eigmin_poisson1D(int *la) {
  double h = 1.0 / ((*la) + 1);
  return 4 * pow(sin(M_PI * h / 2.0), 2);
}

// Calculate optimal alpha for Richardson iteration
double richardson_alpha_opt(int *la) {
  double lambda_min = eigmin_poisson1D(la);
  double lambda_max = eigmax_poisson1D(la);
  return 2.0 / (lambda_min + lambda_max);
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
  double *B = (double *)malloc(*la * sizeof(double));
  double norme_B = cblas_dnrm2(*la, RHS, 1);
  double residu;
  
  *nbite = 0;
  while (*nbite < *maxit) {
      // Copie RHS dans B
      cblas_dcopy(*la, RHS, 1, B, 1);
      // b = b - Ax
      cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, B, 1);
      
      // Calcul du residu
      residu = cblas_dnrm2(*la, B, 1) / norme_B;
      resvec[*nbite] = residu;
      
      // x = x + alpha * (b-Ax)
      cblas_daxpy(*la, *alpha_rich, B, 1, X, 1);
      
      if (residu <= *tol) break;
      (*nbite)++;
  }
  
  printf("\nL'erreur résiduelle obtenu à la fin est de %f\n", residu);
  printf("Le nombre total d'itération pour Richardson alpha est de %d\n", *nbite);
  
  free(B);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {
    
  for (int i = 0; i < *la; i++) {
      for (int j = 0; j < *lab; j++) {
          MB[i * (*lab) + j] = 0.0;
      }
      MB[i * (*lab) + 1] = AB[i * (*lab) + 1];  // Copie diagonale
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {

  for (int i = 0; i < *la; i++) {
      for (int j = 0; j < *lab; j++) {
          MB[i * (*lab) + j] = 0.0;
      }
      MB[i * (*lab) + 1] = AB[i * (*lab) + 1];  // Diagonale
      MB[i * (*lab) + 2] = AB[i * (*lab) + 2];  // Sur-diagonale
  }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite, int *NRHS, int *ipiv, int *info) {
  double *B = (double *)malloc(*la * sizeof(double));
  double norme_B = cblas_dnrm2(*la, RHS, 1);
  double residu;
  int ku_minus = *ku - 1;

  // La factorisation LU pour la matrice bande MB
  dgbtrf_(la, la, kl, &ku_minus, MB, lab, ipiv, info);
  
  *nbite = 0;
  while (*nbite < *maxit) {
      cblas_dcopy(*la, RHS, 1.0, B, 1.0);
      // b - A * x**k
      cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, B, 1);
      
      // Nouvelle erreur résiduelle
      residu = cblas_dnrm2(*la, B, 1) / norme_B;
      resvec[*nbite] = residu;
      
      dgbtrs_("N", la, kl, &ku_minus, NRHS, MB, lab, ipiv, B, la, info, (unsigned long)1);
      cblas_daxpy(*la, 1.0, B, 1, X, 1);

      if (residu <= *tol) break;
      (*nbite)++;
  }
  
  printf("\nL'erreur résiduelle obtenu à la fin est de %f\n", residu);
  printf("Le nombre total d'itération pour Richardson MB est de %d\n", *nbite);
  
    free(B);
}