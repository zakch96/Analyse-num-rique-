#include"iterative_methods.h"

double ** precond_ssor(double ** a, int n, double w){
    //create variables and initialization
    double ** e = create_mat(n);
    double ** e_t = create_mat(n);
    double **  d = create_mat(n);
    double ** d_inv = create_mat(n);
    double ** c = create_mat(n);
    double ** tmp = create_mat(n);
    cout << "in precond_ssor line 11 c = ";
    print_matrix(c, n);


    
    for ( int i = 0; i < n; i++){
        d[i][i] = a[i][i];
        if(a[i][i] != 0)
            d_inv[i][i] = 1/a[i][i];
        for ( int j = 0; j <= i; j++){
            if (j < i){
                e[i][j] = -w *a[i][j];
            }
            else
                e[i][j] = a[i][j];
        }
    }

    e_t = transpose(e, n);
    tmp = dot(d_inv, e_t, n);
    cout << "in precond_ssor line 31 tmp = ";
    print_matrix(tmp, n);
    cout << "in precond_ssor line 33 e_t = ";
    print_matrix(e_t, n);
    c = dot(e, tmp, n);
    cout << "in precond_ssor line 37 c = ";
    print_matrix(c, n);

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            c[i][j] = c[i][j] * 1/(w*(2-w));
            cout << "c[i,j] = " << c[i][j] << endl;
        }
    }



    //free memory
    free_matrix(d,n);
    free_matrix(d_inv,n);
    free_matrix(e,n);
    free_matrix(e_t,n);
    free_matrix(tmp,n);
    // cout << "in precond_ssor c = ";
    // print_matrix(c, n);
    return c;
}

double * grad_conj_precond(double ** a, double * b, int n, double w = 1, double eps = pow(10,-12)){
    //variables initialization
    int max_iter = 7000;
    int k = 0;
    double * x = create_col(n);
    for ( int i = 0; i < n ; i++)
        x[i] = 0.0;
    double * x_1 = create_col(n);
    double * r = create_col(n);
    double * r_1 = create_col(n);
    double * z = create_col(n);
    double * z_1 = create_col(n);
    double * p = create_col(n);
    double * p_1 = create_col(n);
    double * tmp_2 = create_col(n);
    tmp_2 = dot_col(a,x,n);
    double ** c = precond_ssor(a , n, w);
    cout << "precond mat c = ";
    print_matrix(c ,n);
    double alpha = 0.;
    double beta = 0. ;
    //step 0
    for (int i = 0; i < n; i++){
        r[i] = b[i]-tmp_2[i];
    }
    cout << "r = " ;
    print_array(r, n);

    p = solve_cholesky(c, r, n);
    cout << "solve_cholesky solution = ";
    print_array(p,n);

    double * cholesky_test = dot_col(c, p, n);
    cout << "cholesky_test product =" ;
    print_array(cholesky_test, n);

    for (int i = 0; i < n; i++){
        z[i] = p[i];
    }

    for ( int i = 0; i < n; i++){
        x_1[i] = x[i];
        r_1[i] = r[i];
        z_1[i] = z[i];
        p_1[i] = p[i];
    }

    double q = norm(r,n)/norm(b, n);
    cout << "norm(r,n) = " << norm(r,n) << endl;
    cout << "norm(b,n) = " << norm(b,n) << endl;
    cout << "q = " << q << endl;

    //step k
    while( (q >= eps) && k < max_iter ){
        double * tmp = dot_col(a,p,n);
        alpha = scalar_product(r, z, n)/scalar_product(tmp, p, n);
        cout << "alpha = " << alpha  << endl ;
        cout << "x = ";
        print_array(x,n);
        for ( int i = 0; i < n; i++){
            cout << "x_1 = " << x_1[i] << endl ;
            cout << "r_1 = " << x_1[i] << endl ;
            x_1[i] = x_1[i] + alpha * p[i];
            r_1[i] = r_1[i] - alpha * tmp[i];
        }
        z_1  = solve_cholesky(c, r, n);
        double * cholesky_test = dot_col(c, p, n);
        cout << "cholesky_test product =" ;
        print_array(p,n);
        beta = scalar_product(r_1, z_1, n)/scalar_product(r, z, n);
        for(int i = 0; i < n; i++){
            p_1[i] = r_1[i] + beta * p[i];
        }
        for (int i = 0 ; i < n; i++){
            x[i] = x_1[i];
            r[i] = r_1[i];
            z[i] = z_1[i];
            p[i] = p_1[i];
        }
        cout << "r = " ;
        print_array(r, n);
        cout << "b = " ;
        print_array(b, n);
        cout << "norm(r,n) = " << norm(r,n) << endl;
        cout << "norm(b,n) = " << norm(b,n) << endl;
        q = norm(r, n)/norm(b, n);
        cout << "q = " << q << endl; 
        k++;
    }
    cout << "number of iterations k = " << k << endl ; 
    //free memory
    delete [] x_1; 
    delete [] r_1; 
    delete [] r; 
    delete [] p;
    delete [] p_1; 
    delete [] z; 
    delete [] z_1;
    delete [] tmp_2;
    free_matrix(c, n);
    //returning solution
    print_array(x,n);
    cout << "" << endl;
    print_array(x_1,n);

    return x;
}






// double* Grad_Conjug(double** A, double* b, double* x_k, int n)
// {
//   int count_itere = 0;
//   double* r_k = new double[n];
//   double* p_k = new double[n];
//   double* r_k_1 = new double[n];
//   double* x_k_1 = new double[n];
//   double* p_k_1 = new double[n];
//   double* prod_mat_vec_x = dot_col(A, x_k, n);
//   double* prod_mat_vec_p = dot_col(A, p_k, n);
//   for(int i=0; i<n; i++)
//     r_k[i] = b[i] - prod_mat_vec_x[i];
//   for(int i=0; i<n; i++)
//     p_k[i] = r_k[i];
//   double alpha_k;
//   double beta_k;
//   while(1){
//     alpha_k = scalar_product(r_k, r_k, n)/scalar_product(p_k, prod_mat_vec_p, n);
//     for(int i=0; i<n; i++){
//       x_k_1[i] = x_k[i] + alpha_k*p_k[i];
//       r_k_1[i] = r_k[i] - alpha_k*prod_mat_vec_p[i];
//     }
//     count_itere++;
//     if("r_k_1 est suffisamment petit"){
//       cout << "le nombre d'itérations est" << count_itere << endl;
//       return x_k_1;
//     }
//     beta_k = scalar_product(r_k_1, r_k_1, n)/scalar_product(r_k, r_k, n);
//     for(int i=0; i<n; i++)
//       p_k_1[i] = r_k_1[i] + beta_k*p_k[i];
//     for(int i=0; i<n; i++){
//       r_k[i] = r_k_1[i];
//       p_k[i] = p_k_1[i];
//     }
//   }
// }















