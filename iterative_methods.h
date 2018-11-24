#include"cholesky.h"

double ** precond_ssor(double ** a, int n, double w);
double * grad_conj_precond(double ** a, double * b, int n, double w ,double eps );