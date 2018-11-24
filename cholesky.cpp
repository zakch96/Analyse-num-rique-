#include"cholesky.h"


void ldl_decomposition(double **a, int n){
	
	//initialization
	double ** l = create_mat(n);
	double * d = new double [n];
	double * v = new double [n];


	for (int i = 0; i<n; i++){
		d[i] = 0.;
		l[i][i] = 1;
		v[i] = 0;
		for ( int j = 0; j<n; j++){
			if  ( i!=j){
				l[i][j] = 0. ;
			}
		}
	}

	//algorithm
	double sum = 0;
	for ( int i = 0; i<n; i++){
		sum = 0;
		for (int j = 0; j < i-1; j++){
			v[j] = l[i][j]*d[j];
			sum+= l[i][j]*v[j];		
		}
		d[i] = a[i][i] - sum;
		for (int j = i+1; j < n; j++){
			sum = 0;
			for ( int k = 0; k < i-1; k++){
				sum += l[j][k]*v[k];
			}
			l[j][i] = (a[j][i]-sum)/d[i];
		}
	}

	cout << "SHOW L" << endl;
	print_matrix(l,n);
	for (int i = 0; i< n ; i++){
		cout << "d[i] = " << d[i] << endl;
	}
}

void Cholesky_Decomposition(double ** matrix, double ** lower, int n) {

    // Decomposing a matrix into Lower Triangular 
    for (int i = 0; i < n; i++) { 
        for (int j = 0; j <= i; j++) { 
            double sum = 0; 
  
            if (i == j) // summation for diagnols 
            { 
                for (int k = 0; k < j; k++) 
                    sum += pow(lower[j][k], 2); 
                lower[j][j] = sqrt(abs(matrix[j][j] - sum)); 
            } else { 
  
                for (int k = 0; k < j; k++) 
                    sum += (lower[i][k] * lower[j][k]); 
                lower[i][j] = (matrix[i][j] - sum) /lower[j][j]; 
            } 
        } 
    }
    // cout << "Lower Cholesky Matrix = " ;
    // print_matrix(lower,n);
    // cout <<""<<endl;

    // double ** upper = create_mat(n);
    // upper = transpose(lower,n);
    // cout << "" << endl;
    // cout <<"UPPER"<<endl;
    // print_matrix(upper,n);
    // cout << " " << endl;
    // double ** tmp = dot(lower,upper,n);
    // print_matrix(tmp,n);
}

bool element_diagonale_nul(double ** matrix, int n){

	int resultat = false;
	for (int i = 0; i < n; i++){
		for (int j = 0 ; j < n; j++){
			if (matrix[i][j] == 0)
				resultat = true;
		}
	}

  return resultat ;
}


int solve_upper(double ** a, double *x, double * b, int n){
	double sum = 0.;
	int res = 0;

	if (element_diagonale_nul(a, n) == true){
		res = -1;
	}

	x[n-1] = b[n-1]/a[n-1][n-1];
	for (int k = n-2 ; k >=0 ; k--){
		for (int j = k+1; j < n; j++){
			sum+= a[k][j]*x[j];
		}
    x[k] = (b[k] - sum)/a[k][k];
    sum = 0.;
    }

    return res;
}

int solve_lower(double ** a, double *x, double * b, int n){

	double sum = 0.;
	int res = 0;

	if (element_diagonale_nul(a, n) == true){
		res = -1;
	}
	x[0] = b[0]/a[0][0];
	for (int i = 1 ; i < n ; i++){
		for (int j = 0; j < i; j++){
			sum+= a[i][j]*x[j];
		}
		x[i] = (b[i] - sum)/a[i][i];
		sum = 0.;
	}

	return res;
}

double * solve_cholesky(double ** a, double * b, int n){
	int res_solve = 0;
	double * x = new double [n];
	for ( int i = 0 ; i < n ; i++) x[i] = 0. ;
	double ** lower = create_mat(n);
	Cholesky_Decomposition(a, lower,n);
	double ** upper = transpose(lower,n);
	double * y = new double [n];
	for ( int i = 0 ; i < n ; i++) y[i] = 0. ;
	res_solve = solve_lower(lower, y, b, n);
	solve_upper(upper, x, y, n);
	//free_memory
	delete [] y;
	free_matrix(upper, n);
	free_matrix(lower, n);
	return x;
}





