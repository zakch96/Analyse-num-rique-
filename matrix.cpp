#include "matrix.h"

using namespace std;

void print_matrix(double **matrix, int n){
	cout << "[";
	for (int i = 0; i<n ; i++){
		cout << "[" ;
		for (int j = 0; j<n; j++){
			if ( j < n-1)
				cout << matrix[i][j] <<", " ;
			else{
				if (i < n-1)
					cout << matrix[i][j] << "]" << endl ;
				else
					cout << matrix[i][j] << "]]" << endl ;
			}
		}
	}

}

void print_array(double * x, int n){
	cout << "[ ";
	for ( int  i = 0; i < n ; i++ ){
		if ( i <n-1)
			cout << x[i] <<", ";
		else
			cout << x[i] << " ]" << endl;
	}
}

double ** dot(double ** mat1, double ** mat2, int n){
	
	double ** returned_mat = new double *[n];
	double tmp = 0.0;

	for (int i = 0; i <n; i++){
    	returned_mat[i] = new double[n];
    }

    for (int i = 0; i < n; i++){
    	for (int j = 0; j < n; j++){
    		tmp = 0.;
    		for ( int k = 0; k < n; k++){
    			tmp += mat1[i][k]*mat2[k][j];
    		}

    		returned_mat[i][j] = tmp;
    		tmp = 0.0;
    	}
    }

    return returned_mat;
}

double ** dot_diag(double ** mat1, double * diag, int n){

	double ** returned_mat = new double *[n];
	double tmp = 0.0;

	for (int i = 0; i <n; i++){
    	returned_mat[i] = new double[n];
    }

    for (int i = 0; i < n; i++){
    	for (int j = 0; j < n; j++){
    		tmp = mat1[i][j]*diag[j];
    		returned_mat[i][j] = tmp;
   		}
    }

    return returned_mat;
}


double * dot_col(double ** mat, double * col, int n){

	double * returned_col = create_col(n);
	double tmp = 0.0;
	for ( int i = 0; i < n; i++){
		tmp = 0.0;
		for (int k = 0; k < n ; k++){
			tmp+= mat[i][k] * col[k];
		}
		returned_col[i] = tmp;
	}
	return returned_col;
}

double ** sum_mat(double ** mat1, double ** mat2, int n){

	double ** returned_mat = new double *[n];
	double tmp = 0.0;

	for (int i = 0; i <n; i++){
    	returned_mat[i] = new double[n];
    }
    for ( int i = 0; i < n; i++){
    	for ( int j = 0; j < n; j++){
    		tmp = mat1[i][j] + mat2[i][j];
    		returned_mat[i][j] = tmp;
    	}
    }

    return returned_mat;
}

double ** create_mat(int n){
	
	double ** returned_mat = new double *[n];
	
	for (int i = 0; i <n; i++){
    	returned_mat[i] = new double[n];
    }

    for ( int k = 0; k < n ; k++){
    	for (int l = 0; l < n; l++){
    		returned_mat[k][l] = 0.0;
    	}
    }

    return returned_mat;
}

double * create_col(int n){
	double * col = new double [n];
	for (int i = 0; i < n; i++)
		col[i] = 0. ;
	return col;
}

void mat_initialize(double ** mat, int n, int m, double a, double b){
	int dim = n*m;
	double h = a/m+1;
	double k = b/n+1;
	double d = (2./pow(h,2)) + (2./pow(k,2)) ; 

	for (int i = 0; i < dim ; i++ ){
		if (i >= 1)
			mat[i][i-1] = -1./pow(h,2); 
		if (i <= dim-2)
			mat[i][i+1] = -1/pow(h,2) ;
		if (i <= dim-2)
			mat[i][i+1] = -1./pow(h,2);
		if (i >= n)
			mat[i][i-n] = -1./pow(k,2) ;
		if ( i <= dim - n - 1)
			mat[i][i+n] = -1./pow(k,2);
		mat[i][i] = d ;
	}

	for (int j = m; j < dim-m+1; j=j+m){
		mat[j][j-1] = 0. ;
	}
	for ( int k = dim - m-1; k>m-2; k = k-m){
		mat[k][k+1] = 0. ;
	} 
}

void column_initialize(double * col, int n, int m, double a, double b){
	//VARIABLES INITIALIZATION
	//int dim = n * m;
	//int ind = 0;
	double * x = new double[m];
	double * y = new double[n];
	double h = a/m+1;
	double k = b/n+1;

	for(int i = 1; i < m; i++){
		for (int j = 0; j < n; j++)
			col[i+(j-1)*m] = f(i*h, j*h);
	}

	for ( int i = 1; i < m-1; i++){
		col[i + (n-1)*m] = f(i, n) + (g(x[i], b)/pow(k,2));
	}
	for ( int i = 1; i < m-1; i++){
		col[i] = f(i, 1) + (g(x[i], 0)/pow(k,2));
	}
	for ( int j = 1; j < n-1; j++){
		col[1+(j-1)*m] = f(1, j) + (g(0, y[j])/pow(h,2));
	}
	for ( int j = 1; j < n-1; j++){
		col[m+(j-1)*m] = f(m, j) + (g(a, y[j])/pow(h,2));
	}
	col[0] = f(h, k) + g(h, 0)/pow(k, 2) + g(0, k)/pow(h, 2);
	col[m-1] = f(h*m, k) + g(m*h, 0)/pow(k, 2) + g(a, k)/pow(h, 2);
	col[(n-1)*m-1] = f(h, n*k) + g(0, n*k)/pow(h, 2) + g(h, b)/pow(k, 2);
	col[m + (n-1)*m-1] = f(m*h, n*k) + g(a, n*k)/pow(h, 2) + g(m*h, b)/pow(k, 2);
	return;
}

double scalar_product(double *col1, double * col2, int n){
	double res = 0.0 ; 
	for ( int i = 0; i < n; i++){
		res += col1[i] * col2[i];
	}
	return res;
}

double norm(double * col1, int n){
	double res = 0.0;
	for (int i = 0; i < n; i++){
		res += col1[i]*col1[i];
	}
	res = sqrt(res);
	return res;
}

double ** transpose(double ** mat, int n){
	double ** transpose_mat = create_mat(n);
	for(int i = 0; i < n; i++){
		for ( int j = 0; j < n; j++)
			transpose_mat[i][j] = mat[j][i];
	}

	return transpose_mat;
}


void free_matrix(double ** mat, int n){
	for(int i = 0; i < n; ++i)
    	delete [] mat[i];
	delete [] mat;
}


