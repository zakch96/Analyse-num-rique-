#include"func.h"
#include"iterative_methods.h"
#include"eigen_value_computation.h"

int main(){
	double **mat1 = create_mat(5) ;
	double **mat2 = create_mat(5) ;
	double **mat3 = create_mat(5) ;
	double *col = new double[5];

	// mat_initialize(mat1, 5, 5, 2, 5);
	// print_matrix(mat1, 9);

	// dot_col(mat1,col,5);
	// sum_mat(mat1,mat2,5);
	// dot(mat1,mat2,5);

	
    
    delete[] col;

    double a [5][5] = {{320, 4, 5, 6, 10},
    				   {4, 170, 7 ,71,20},
    				   {5, 7, 200, 5,  2}, 
    				   {6, 71,5,350,4}, 
    				   {10, 20,2,4,190}};

    double a2 [4][4] = {{1, 1, 1, 1 },
                        {1, 5, 5, 5 },
                        {1, 5, 14,14}, 
                        {1, 5,14,15}};
    double a3[4][4] = {{1, 2, 3, 4},
                       {0, 3, 5, 7},
                       {0, 0, 7, 9},
                       {0, 0, 0, 4}};

    double ** lower1 = create_mat(5);
    double ** lower2 = create_mat(4);

    
    for (int i = 0 ; i < 5; i++){
    	for (int j = 0; j < 5; j++){
    		mat1[i][j] = a[i][j] ;
    	}
    }
    for (int i = 0 ; i < 4; i++){
        for (int j = 0; j < 4; j++){
            mat2[i][j] = a2[i][j] ;
            mat3[i][j] = a3[i][j];
        }
    }
    // cout << "print before" << endl;
    // print_matrix(mat1, 5);
    // ldl_decomposition(mat1 , 5);
    // cout << "after calling solve_cholesky" << endl;
    // print_matrix(mat1, 5);
    // Cholesky_Decomposition(mat1, lower1, 5);
    // Cholesky_Decomposition(mat2, lower2, 4);
    double arr [4] = {1,2,3,3};
    double * b = new double[4];
    for ( int i = 0 ; i < 5 ; i++) b[i] = arr[i];
    cout << "b = " ;
    print_array(b, 4);
    double * sol = solve_cholesky(mat2, b, 4);
    cout << "sol = " ;
    print_array(sol, 4);
    double * solution = create_col(4);
    solution = dot_col(mat2, sol, 4);
    cout << "last print solution b = " ;
    print_array(solution, 4);
    //double w = 1.5;
    double ** c = precond_ssor(mat2, 4, 1);
    cout << "PRECOND SSOR TEST" << endl;
    print_matrix(c, 4);
    //grad conj test
    cout << "Grad conj SSOR test " << endl;
    // double * b2 = create_col(4);
    // b2 = grad_conj_precond(mat2, b, 4, 0.5, pow(10.,-12));
    // double * tmp2= dot_col(mat2, b2,4);
    // cout << "tmp2 = " ; 
    // print_array(tmp2, 4);
    // cout << "b2 =";
    // print_array(b2, 4);
    cout << "" << endl ; 
    cout << "MAX EIGEN VALUE TEST " << endl;
    double * x = create_col(4);
    double * x2 = create_col(4);
    double x_copy[] = {1.5,-1,0,1};
    for (int i = 0; i < 4; i++) x[i] = x_copy[i];
    double max_lambda = max_eigenvalue(mat3, 4, x, 0.0001);
    cout << "max_lambda = " << max_lambda << endl;
    cout << "MIN EIGEN VALUE TEST " << endl;
    double x_copy2[] = {1.5,-1,0,1};
    for (int i = 0; i < 4; i++) x2[i] = x_copy2[i];

    double min_lambda = min_eigenvalue(mat3, 4, x2, 0.0001);
    cout << "min_lambda = " << min_lambda << endl;

    //free memory
    free_matrix(mat1,5);
    free_matrix(mat2,4);
    free_matrix(mat3,4);
    free_matrix(lower1,5);
    free_matrix(lower2,4);

	return 0;
}