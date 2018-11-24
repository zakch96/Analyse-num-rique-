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