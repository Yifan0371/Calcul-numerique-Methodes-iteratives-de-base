/******************************************/
/* tp2_poisson1D_direct.c                 */
/* Main program to solve the 1D Poisson   */
/* problem using LU factorization         */
/******************************************/
#include "lib_poisson1D.h"
#include "time.h"

#define TRF 0  // Option for LAPACK DGBTRF
#define TRI 1  // Option for custom LU decomposition (tridiagonal)
#define SV  2  // Option for LAPACK DGBSV

int main(int argc, char *argv[])
{
  struct timespec start, end; // For time measurement
  double elapsed_time;        // Time storage variable
  int ierr, info = 1;


  int nbpoints, la;           // Number of points and discretized size
  int ku, kl, kv, lab;        // Matrix band parameters
  int *ipiv;                  // Pivot indices for LU factorization
  int NRHS = 1;               // Number of RHS vectors
  int IMPLEM = 0;             // Implementation choice




  double T0 = -5.0, T1 = 5.0; // Boundary conditions
  double *RHS, *EX_SOL, *X;   // RHS, exact solution, grid points
  double *AB;                 // Band matrix
  double relres;              // Relative forward error



  // Handle optional argument to choose the implementation method
  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  /* Initialization */
  nbpoints = 10;
  la = nbpoints - 2; // Internal grid points excluding boundaries
  RHS = (double *)malloc(sizeof(double) * la);
  EX_SOL = (double *)malloc(sizeof(double) * la);
  X = (double *)malloc(sizeof(double) * la);

  printf("--------- Poisson 1D ---------\n\n");

  // Generate grid points, RHS vector, and exact solution
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

  // Write vectors to files for verification
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  /* Matrix band parameters */
  kv = 1;                     // Main diagonal position
  ku = 1;                     // Number of super-diagonals
  kl = 1;                     // Number of sub-diagonals
  lab = kv + kl + ku + 1;     // Leading dimension for band storage

  AB = (double *)malloc(sizeof(double) * lab * la);
  ipiv = (int *)calloc(la, sizeof(int));

  // Set up the band matrix for Poisson 1D
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  /* Solve using different LU implementations */
  if (IMPLEM == TRF) {
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info); // LAPACK DGBTRF
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Execution time for DGBTRF: %f seconds\n", elapsed_time);
  }

  if (IMPLEM == TRI) {
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info); // Custom LU for tridiagonal matrix
  }

  if (IMPLEM == TRI || IMPLEM == TRF) {
    if (info == 0) {
      clock_gettime(CLOCK_MONOTONIC_RAW, &start);
      dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info); // Solve triangular system
      clock_gettime(CLOCK_MONOTONIC_RAW, &end);
      elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
      printf("Execution time for DGBTRS: %f seconds\n", elapsed_time);
    }
  }

  if (IMPLEM == SV) {
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info); // Solve using LAPACK DGBSV
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Execution time for DGBSV: %f seconds\n", elapsed_time);
  }

  // Write the LU factorized matrix and solution to files
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
  write_xy(RHS, X, &la, "SOL.dat");

  /* Compute and display relative forward error */
  relres = relative_forward_error(RHS, EX_SOL, &la);
  printf("\nThe relative forward error is relres = %e\n", relres);

  /* Free memory */
  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);

  printf("\n\n--------- End -----------\n");
  return 0;
}
