double* Grad_Conjug_Precond(double** A, double* b, double* x_k, int n, double epsilon, double w)
{
  int count_itere = 0;
  double* r_k = new double[n];
  double* p_k = new double[n];
  double* z_k = new double[n];
  double* r_k_1 = new double[n];
  double* x_k_1 = new double[n];
  double* p_k_1 = new double[n];
  double** C = precond_ssor(A, n, w);
  double* prod_mat_vec_x = dot_col(A, x_k, n);
  for(int i=0; i<n; i++)
    r_k[i] = b[i] - prod_mat_vec_x[i];
  double* prod_mat_vec_p_0 = solve_cholesky(C, r_k, n);
  for(int i=0; i<n; i++){
    p_k[i] = prod_mat_vec_p_0[i];
    z_k[i] = prod_mat_vec_p_0[i];
  }
  double* prod_mat_vec_p = new double[n];
  double alpha_k;
  double beta_k;
  while(1){
    prod_mat_vec_p = dot_col(A, p_k, n);
    alpha_k = scalar_product(r_k, r_k, n)/scalar_product(p_k, prod_mat_vec_p, n);
    for(int i=0; i<n; i++){
      x_k_1[i] = x_k[i] + alpha_k*p_k[i];
      r_k_1[i] = r_k[i] - alpha_k*prod_mat_vec_p[i];
    }
    count_itere++;
    if(scalar_product(r_k_1, r_k_1, n) < epsilon){
      delete r_k;
      delete p_k;
      delete r_k_1;
      delete p_k_1;
      delete prod_mat_vec_x;
      delete prod_mat_vec_p;
      cout << "le nombre d'itérations est :   " << count_itere << endl;
      return x_k_1;
    }
    prod_mat_vec_p_0 = solve_cholesky(C, r_k, , n);
    for(int i=0; i<n; i++)
      z_k[i] = prod_mat_vec_p_0[i];
    beta_k = scalar_product(r_k_1, r_k_1, n)/scalar_product(r_k, r_k, n);
    for(int i=0; i<n; i++)
      p_k_1[i] = r_k_1[i] + beta_k*p_k[i];
    for(int i=0; i<n; i++){
      r_k[i] = r_k_1[i];
      p_k[i] = p_k_1[i];
      x_k[i] = x_k_1[i];
    }
  }
}