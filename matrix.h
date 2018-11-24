#include <iostream>
#include<cmath>
#include<vector>
#include<string>

using namespace std;

void print_matrix(double **matrix, int n);

void print_array(double * x, int n);

double ** dot(double ** mat1, double ** mat2, int n);

double * dot_col(double ** mat, double * col, int n);

double ** sum_mat(double ** mat1, double ** mat2, int n);

double ** create_mat(int n);

double * create_col(int n);

void mat_initialize(double ** mat, int n, int m, double a, double b);

double ** transpose(double ** mat, int n);

void free_matrix(double ** mat, int n);

double scalar_product(double *col1, double * col2, int n);

double norm(double * col1, int n);

// void column_initialize(double * col, int m, int n, double (&f)(double x, double y), double (&g)(double k, double l));

double ** dot_diag(double ** mat1, double * diag, int n);
