#include"matrix.h"

void ldl_decomposition(double **a, int n);
void Cholesky_Decomposition(double ** matrix, double ** lower, int n);
int solve_upper(double ** a, double *x, double * b, int n);
int solve_lower(double ** a, double *x, double * b, int n);
bool element_diagonale_nul(double ** matrix, int n);
double * solve_cholesky(double ** a, double * b, int n);