/**********************************************/
/* lib_poisson1D.c                            */
/* Bibliothèque numérique pour résoudre       */
/* l'équation de Poisson 1D (chaleur)         */
/**********************************************/
#include "lib_poisson1D.h"
#include <cblas.h> /* pour cblas_daxpy, cblas_dnrm2 */

/* Initialiser la matrice de Poisson au format bande (colonne-major) */
void set_GB_operator_colMajor_poisson1D(double *AB, int *lab, int *la, int *kv)
{
    int n   = *la;  /* Nombre de points intérieurs */
    int ld  = *lab; /* Nombre de lignes dans la matrice bande */
    int k   = *kv;  /* Position du pivot principal */

    /* Remplir la matrice avec des zéros */
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < ld; i++) {
            AB[i + j*ld] = 0.0;
        }
    }

    /* Construire la matrice tridiagonale pour l'équation de Poisson */
    for (int j = 0; j < n; j++) {
        /* Diagonale principale */
        AB[k+1 + j*ld] = 2.0;
        /* Sous-diagonale (commence à la deuxième colonne) */
        if (j > 0) {
            AB[k + j*ld] = -1.0;
        }
        /* Sur-diagonale (jusqu'à l'avant-dernière colonne) */
        if (j < n-1) {
            AB[k+2 + j*ld] = -1.0;
        }
    }
}

/* Initialiser le second membre (RHS) avec des conditions de Dirichlet */
void set_dense_RHS_DBC_1D(double *RHS, int *la, double *BC0, double *BC1)
{
    int n = *la; /* Nombre de points intérieurs */
    RHS[0]    = *BC0;
    for (int i = 1; i < n - 1; ++i) {
        RHS[i] = 0.0;
    }
    RHS[n - 1] = *BC1;
}

/* Calculer la solution analytique pour le problème de Poisson 1D */
void set_analytical_solution_DBC_1D(double *EX_SOL, double *X, int *la, double *BC0, double *BC1)
{
    int n = *la; /* Nombre de points intérieurs */
    double DT = (*BC1 - *BC0);
    for (int i = 0; i < n; i++) {
        EX_SOL[i] = *BC0 + X[i] * DT;
    }
}

/* Générer un maillage uniforme pour le domaine */
void set_grid_points_1D(double *x, int *la)
{
    int n = *la; /* Nombre de points intérieurs */
    double h = 1.0 / (n + 1);
    for (int jj = 0; jj < n; jj++) {
        x[jj] = (jj + 1) * h;
    }
}

/* Calculer l'erreur relative entre deux vecteurs */
double relative_forward_error(double *x, double *y, int *la)
{
    int n = *la; /* Nombre de points intérieurs */
    cblas_daxpy(n, -1.0, x, 1, y, 1); /* Soustraction des vecteurs */
    return cblas_dnrm2(n, y, 1) / cblas_dnrm2(n, x, 1); /* Calcul de l'erreur relative */
}

/* Factorisation LU pour les matrices tridiagonales au format bande */
/* Hypothèse: la matrice AB est une matrice tridiagonale */
int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info)
{
    int N   = *la;  /* Nombre de points intérieurs */
    int ld  = *lab; /* Nombre de lignes dans la matrice bande */
    int i_kv = (*ku); /* Position du pivot principal */

    /* Factorisation LU colonne par colonne */
    for (int j = 1; j < N; j++) {
        AB[i_kv + j*ld] = AB[i_kv + j*ld] / AB[i_kv+1 + (j-1)*ld]; /* L_j */
        AB[i_kv+1 + j*ld] -= AB[i_kv + j*ld]*AB[i_kv+2 + (j-1)*ld]; /* Mise à jour de la diagonale principale */
    }

    /* Initialiser le tableau des permutations */
    for (int ip = 0; ip < N; ip++) {
        ipiv[ip] = ip + 1;
    }
    *info = 0;
    return 0;
}
