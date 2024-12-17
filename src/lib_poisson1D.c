/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library for solving 1D Poisson   */
/* equation (Heat equation)                   */
/**********************************************/
#include "lib_poisson1D.h"

/* Initialize the Poisson matrix in band format (column-major) */
void set_GB_operator_colMajor_poisson1D(double *AB, int *lab, int *la, int *kv)
{
  int rl = *la, cl = *lab, col0 = *kv, i = 1;

  while (i < rl * cl - 1)
  {
    if (i % cl == 0) { for (int j = 0; j < col0; ++j) AB[i++] = 0; }
    if (i == col0) AB[i] = 0;
    else if ((i % cl == col0) || (i % cl == col0 + 2)) AB[i] = -1;
    else if (i % cl == col0 + 1) AB[i] = 2;
    i++;
  }
  AB[rl * cl - 1] = 0;
}

/* Set the right-hand side (RHS) vector with Dirichlet boundary conditions */
void set_dense_RHS_DBC_1D(double *RHS, int *la, double *BC0, double *BC1)
{
  RHS[0] = *BC0;
  for (int i = 1; i < *la - 1; ++i) RHS[i] = 0;
  RHS[*la - 1] = *BC1;
}

/* Compute the analytical solution for the 1D Poisson problem */
void set_analytical_solution_DBC_1D(double *EX_SOL, double *X, int *la, double *BC0, double *BC1)
{
  const double DT = *BC1 - *BC0;
  for (int i = 0; i < (*la); i++) EX_SOL[i] = *BC0 + X[i] * DT;
}

/* Generate a uniform grid for the domain */
void set_grid_points_1D(double *x, int *la)
{
  double h = 1.0 / ((*la) + 1);
  for (int jj = 0; jj < (*la); jj++) x[jj] = (jj + 1) * h;
}



/* Compute the relative forward error between two vectors */
double relative_forward_error(double *x, double *y, int *la)
{
  cblas_daxpy(*la, -1, x, 1, y, 1);
  return cblas_dnrm2(*la, y, 1) / cblas_dnrm2(*la, x, 1);
}



/* LU factorization for tridiagonal matrices in band format */
int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info)
{
  for (int ii = *lab - 1; ii < *lab * *la; ii++) {
    if (!(ii % *lab)) {
      if (AB[ii - 2]) AB[ii - 1] = AB[ii + 1] / AB[ii - 2];
      AB[ii + 2] -= AB[ii - 1] * AB[ii + 1];
    }
  }

  for (int ip = 0; ip < *la; ip++) ipiv[ip] = ip + 1;
  *info = 0;
  return 0;
}