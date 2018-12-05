#include"func.h"
#include"iterative_methods.h"
#include"eigen_value_computation.h"

int main(){
	double **mat1 = create_mat(5) ;
	double **mat2 = create_mat(5) ;
	double **mat3 = create_mat(5) ;

	

    double a [5][5] = {{320, 4, 5, 6, 10},
    				   {4, 300, 7 ,71,20},
    				   {5, 7, 200, 5,  2}, 
    				   {6, 71,5,350,4}, 
    				   {10, 20,2,4,250}};

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
            mat2[i][j] = a2[i][j];
            mat3[i][j] = a3[i][j];
        }
    }
    // print_matrix(mat1, 5);
    // ldl_decomposition(mat1 , 5);
    // cout << "after calling solve_cholesky" << endl;
    // print_matrix(mat1, 5);
    // Cholesky_Decomposition(mat1, lower1, 5);
    // double * col_1 = create_col(5);
    // double col1[5] =  {3, 8, 1, 2, 7};
    // for ( int i = 0; i < 5 ; i++) col_1[i] = col1[i];
    // double * sol_1 = solve_cholesky(mat1, col_1, 5);
    // cout << "TEST CHOLESKY CAS = 5" << endl;
    // cout << "cholesky_sol = " ;
    // double *product_1 = create_col(5);
    // product_1 = dot_col(mat1, sol_1, 5);
    // print_array(product_1, 5);
 

    // double arr [4] = {1,2,3,3};
    // double * b = new double[4];
    // for ( int i = 0 ; i < 5 ; i++) b[i] = arr[i];
    // cout << "b = " ;
    // print_array(b, 4);
    // double * sol = solve_cholesky(mat2, b, 4);
    // cout << "sol = " ;
    // print_array(sol, 4);
    // double * solution = create_col(4);
    // solution = dot_col(mat2, sol, 4);
    // cout << "last print solution b = " ;
    // print_array(solution, 4);
    // //double w = 1.5;
    // cout << "PRECOND SSOR TEST" << endl;
    // double ** c = precond_ssor(mat2, 4, 0.5);
    // print_matrix(c, 4);
    // double ** c_anss = Create_Precond_C(mat2, 4, 0.5);
    // cout << "anass _ matrix = " << endl;
    // print_matrix(c_anss, 4);
    // //grad conj test
    // cout << "Grad conj SSOR test " << endl;
    // // double * b2 = create_col(4);
    // // b2 = grad_conj_precond(mat2, b, 4, 0.5, pow(10.,-12));
    // // double * tmp2= dot_col(mat2, b2,4);
    // // cout << "tmp2 = " ; 
    // // print_array(tmp2, 4);
    // // cout << "b2 =";
    // // print_array(b2, 4);
    // cout << "" << endl ; 
    // cout << "MAX EIGEN VALUE TEST " << endl;
    // double * x = create_col(4);
    // double * x2 = create_col(4);
    // double x_copy[] = {1.5,-1,0,1};
    // for (int i = 0; i < 4; i++) x[i] = x_copy[i];
    // double max_lambda = max_eigenvalue(mat3, 4, x, 0.0001);
    // cout << "max_lambda = " << max_lambda << endl;
    // cout << "MIN EIGEN VALUE TEST " << endl;
    // double x_copy2[] = {1.5,-1,0,1};
    // for (int i = 0; i < 4; i++) x2[i] = x_copy2[i];

    // double min_lambda = min_eigenvalue(mat3, 4, x2, 0.0001);
    // cout << "min_lambda = " << min_lambda << endl;

    //SYSTEM RESOLUTION
    cout << "" << endl;
    cout << "SYSTEM RESOLUTION TEST" << endl;
    cout << "" << endl;

    int n = 6;
    int m = 5;
    int dim = n*m;
    double a0 = 1.0;
    double b0 = 2.0;
    double ** mat = create_mat(dim);
    mat_initialize(mat, n, m, a0, b0);
    double symetrique = 0.0;
    for (int i = 0; i < dim ; i++){
        for(int j = 0; j < dim; j++){
            if (mat[i][j] != mat[j][i]){
                cout <<"INITIALIZED MATRIX IS NOT SYMMETRIC" << endl;
                symetrique = -1;
            }
        }
    }
    cout << "symetrique = " <<  symetrique << endl;
    double * col = create_col(dim);
    column_initialize(col, n, m, a0, b0);
    cout << " " << endl;
    cout << "after initialization__ b = ";
    print_array(col, dim);
    double w = 0.5;
    double eps = 0.01;
    //double ** c2 = precond_ssor(mat, dim, w);
    double * b2 = grad_conj_precond(mat, col, dim, w, eps);
    cout << " " << endl;
    cout << "after calling grad conj,  solution b2 = ";
    print_array(b2, dim);
    cout << "col = " ;
    print_array(col, dim);
    cout << "prod=" ;
    double * tmp__0 = dot_col(mat, b2, dim);
    print_array(tmp__0,dim);

    cout << "TESTING WITH CHOLESKY " << endl;
    double * just_testing = solve_cholesky(mat,col , dim);
    cout << "test = ";
    print_array(just_testing, dim);
    double * cholesky_solution = dot_col(mat, just_testing, dim);
    cout << "(must equal) dot_solution = " ;
    print_array(cholesky_solution, dim);
    cout << " " << endl;
    cout << "after initialization_col = ";
    print_array(col, dim);
    cout << " " << endl;
    cout << "after initialization_mat = ";
    print_matrix(mat, dim);

    // cout << "GRAD CONJ TEST " << endl;
    // double * x_k= create_col(dim);
    // double * grad_conj_test = Grad_Conjug(mat, col, x_k, dim, 0.001);
    // double * grad_conj_dot = dot_col(mat, grad_conj_test, dim);
    // cout << "grad_conj_dot = " ;
    // print_array(grad_conj_dot, n);

    //TESTING GRAD COJUG
    cout << "testing grad conjug" << endl;
    cout << " " << endl;

    double * x_conjug  = create_col(dim);
    double * grad_conjug_test = Grad_Conjug(mat, col, x_conjug, dim,eps);
    double * dot_conjug = dot_col(mat, x_conjug, dim);
    cout << "dot_conjug = ";
    print_array(dot_conjug, dim);
    double err = erreur_2(grad_conjug_test, n, m, a0, b0);
    cout << "err_grad = " << err  << endl;
    double err_cholesky = erreur_2(just_testing, n, m, a0, b0);
    cout << "err_cholesky = " << err_cholesky  << endl;

    


    //free memory
    //delete[] col_1;
    free_matrix(mat1,5);
    free_matrix(mat2,4);
    free_matrix(mat3,4);
    free_matrix(lower1,5);
    free_matrix(lower2,4);




	return 0;

}