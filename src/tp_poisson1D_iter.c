/******************************************/
/* tp2_poisson1D_iter.c                   */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"

#define ALPHA 0
#define JAC 1
#define GS 2

int main(int argc, char *argv[])
{
    int ierr;
    int jj;
    int nbpoints, la;
    int ku, kl, lab, kv;
    int *ipiv = NULL;  
    int info;
    int NRHS;
    int IMPLEM = 0;
    double T0, T1;
    double *RHS, *SOL, *EX_SOL, *X;
    double *AB;
    double *MB;
    
    double temp, relres;
    double opt_alpha;

    
    if (argc == 2) {
        IMPLEM = atoi(argv[1]);
    } else if (argc > 2) {
        perror("Application takes at most one argument");
        exit(1);
    }
   

    /* Size of the problem */
    NRHS = 1;
    nbpoints = 12;
    la = nbpoints - 2;

    /* Dirichlet Boundary conditions */
    T0 = 5.0;
    T1 = 20.0;

    printf("--------- Poisson 1D ---------\n\n");

  
    RHS    = (double *) malloc(sizeof(double) * la);
    SOL    = (double *) calloc(la, sizeof(double));
    EX_SOL = (double *) malloc(sizeof(double) * la);
    X      = (double *) malloc(sizeof(double) * la);

    /* Setup the Poisson 1D problem */
    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    
    write_vec(RHS,    &la, "RHS.dat");
    write_vec(EX_SOL, &la, "EX_SOL.dat");
    write_vec(X,      &la, "X_grid.dat");

    kv = 0;
    ku = 1;
    kl = 1;
    lab = kv + kl + ku + 1;

    AB = (double *) malloc(sizeof(double) * lab * la);
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

   
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

    /********************************************/
    /* 1. Solution (Richardson with optimal alpha) */

    /* Computation of optimum alpha */
    opt_alpha = richardson_alpha_opt(&la);
   
    printf("Optimal alpha for simple Richardson iteration is : %lf\n", opt_alpha);

    double tol   = 1e-3;
    int    maxit = 1000;
    double *resvec;
    int    nbite = 0;

    resvec = (double *) calloc(maxit, sizeof(double));

    if (IMPLEM == ALPHA) {
        richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
      
        write_vec(resvec, &nbite, "RESVEC_alpha.dat");
    }

    /********************************************/
    /* 2. Richardson General Tridiag + Precondition */

    kv = 1;  
    ku = 1;
    kl = 1;
    MB = (double *) malloc(sizeof(double) * lab * la);

    
    ipiv = (int *) malloc(la * sizeof(int));

    if (IMPLEM == JAC) {
        extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
    } else if (IMPLEM == GS) {
        extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
    }

    if (IMPLEM == JAC || IMPLEM == GS) {
        write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "MB.dat");
        richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit,
                      resvec, &nbite, &NRHS, ipiv, &info);
       
        if (IMPLEM == JAC) {
            write_vec(resvec, &nbite, "RESVEC_jacobi.dat");
        } else {
            write_vec(resvec, &nbite, "RESVEC_GS.dat");
        }
    }

 
    write_vec(SOL, &la, "SOL.dat");

   
    write_vec(resvec, &nbite, "RESVEC.dat");

    free(resvec);
    free(RHS);
    free(SOL);
    free(EX_SOL);
    free(X);
    free(AB);
    free(MB);
    if (ipiv) free(ipiv);

    printf("\n\n--------- End -----------\n");
    return 0;
}
