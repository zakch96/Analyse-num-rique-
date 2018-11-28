#include"eigen_value_computation.h"

double max_eigenvalue(double ** a, int n, double * x, double eps){
	//initialization
	double old_val = 0.0;
	double new_val = 1.0;
	double * x_1 = create_col(n);
	double * y = create_col(n);
	int k = 0;

	while (abs(new_val - old_val) >= eps) {
		old_val = new_val;
		print_matrix(a , n);
		print_array(x,n);
		y = dot_col(a, x, n);
		print_array(y, n);
		for ( int i = 0; i < n ; i++){
			x_1[i] = y[i]/norm(y, n);
		}
		new_val = scalar_product(x, y, n);
		k++;
		for ( int i = 0; i < n; i++)
			x[i] = x_1[i];
	}
	cout << "number of iterations = " << k << endl;

	// free memory
	delete [] x_1;
	delete [] y;

	return new_val;
}

double min_eigenvalue(double ** a, int n, double * x, double eps){
	//initialization
	double old_val = 0.5;
	double new_val = 1.0;
	double * u = create_col(n);
	//double * y = create_col(n);
	double * x_1 = create_col(n);
	int k = 0;
	print_matrix(a , n);
	print_array(x,n);

	while (abs(new_val - old_val) >= eps ){
		cout << "new_val:  " << new_val << endl;
		cout << "old_val:  " << old_val << endl;

		print_matrix(a , n);
		print_array(x,n);
		old_val = new_val;
		u = solve_cholesky(a, x, n);
		for(int i = 0; i < n; i++)
			x_1[i] = u[i]/norm(u, n);
		new_val = 1/scalar_product(x, u, n);
		for ( int i = 0; i < n; i++)
			x[i] = x_1[i];
	}
	cout << "number of iterations = " << k << endl;

	delete [] u;
	//delete [] y;
	delete [] x_1;

	return new_val;
}


