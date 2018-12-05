#include"cholesky.h"

double ** precond_ssor(double ** a, int n, double w);
double * grad_conj_precond(double ** a, double * b, int n, double w ,double eps );
double** Create_Precond_C(double** A, int n, double w);
double* Grad_Conjug(double** A, double* b, double* x_k, int dim, double epsilon) ;
double* Grad_Conjug_Precond(double** A, double* b, double* x_k, int n, double epsilon, double w);
double erreur_2(double * x, int n, int m, double a, double b);