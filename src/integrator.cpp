#include "integrator.hpp"


int dy_general_model(double t, const double y[], double dy[], void* params) {
	int status = GSL_SUCCESS; // return value

	OdeParamStruct* pars = (OdeParamStruct*) params;
	double* C = pars->C; // contact matrix
	double beta = pars->beta;
	double gamma = pars->gamma;
	double kappa = pars->kappa;
	int na = pars->na;
	int ne = pars->ne;
	int ni = pars->ni;
	
	/* NB: matrix views of arrays work as follows:
	 * C_mat(i,j) = C[i*n2 + j], where n2 is the number of rows.
	 * Suppose that C = [1,2,3,4,5,6], and n2 = 3, then
	 * C_mat = [1,2,3;4,5,6]
	 */
	gsl_matrix_const_view C_mat = gsl_matrix_const_view_array(C, na, na);
	
	// make y vector views
	gsl_vector_const_view S_vec = gsl_vector_const_view_array(y, na);
	/* gsl_matrix_const_view_array has arguments: array, # rows, # columns
	 * E_vecs has as rows E_1, E_2, ..., E_ne (if ne > 0)
	 * I_vecs has as rows I_1, I_2, ..., E_ni
	 */
	gsl_matrix_const_view I_vecs = gsl_matrix_const_view_array(y+(1+ne)*na, ni, na);
	
	// make dy vector views
	gsl_vector_view dS_vec = gsl_vector_view_array(dy, na);
	gsl_matrix_view dI_vecs = gsl_matrix_view_array(dy+(1+ne)*na, ni, na);
	
	// make the dy vector zero
	gsl_vector_set_zero(&dS_vec.vector);
	gsl_matrix_set_zero(&dI_vecs.matrix);
		
	// dS/dt = -beta diag(S) C (I_1 + I_2 + ... + I_ni)
	gsl_vector* I_sum = gsl_vector_calloc(na); // initiated zero
	for ( int i = 0; i < ni; ++i ) {
		gsl_vector_const_view I = gsl_matrix_const_row(&I_vecs.matrix, i);
		status = gsl_blas_daxpy(1.0, &I.vector, I_sum);		
		if ( status != GSL_SUCCESS ) return status;
	}
	status = gsl_blas_dgemv(CblasNoTrans, -beta, &C_mat.matrix, I_sum, 0.0, &dS_vec.vector);
	if ( status != GSL_SUCCESS ) return status;
	status = gsl_vector_mul(&dS_vec.vector, &S_vec.vector);
	if ( status != GSL_SUCCESS ) return status;
	
	if ( ne > 0 ) { // there is at least one E-compartment
		// repeat some steps taken for the I compartment (but then for E)
		gsl_matrix_const_view E_vecs = gsl_matrix_const_view_array(y+na, ne, na);
		gsl_matrix_view dE_vecs = gsl_matrix_view_array(dy+na, ne, na);
		gsl_matrix_set_zero(&dE_vecs.matrix);

		/* dE_1/dt = -dS/dt - ne*kappa E_1
		 * dE_2/dt = ne*kappa E_1 - ne*kappa E_2
		 * ...
		 */
		for ( int i = 0; i < ne; ++i ) {
			if ( i == 0 ) {
				// select rows from the dE and E matrices
				gsl_vector_view dE0 = gsl_matrix_row(&dE_vecs.matrix, 0);
				gsl_vector_const_view E0 = gsl_matrix_const_row(&E_vecs.matrix, 0);
				// compute dE
				status = gsl_blas_daxpy(-1.0, &dS_vec.vector, &dE0.vector); // y -> y + a x
				if ( status != GSL_SUCCESS ) return status;
				status = gsl_blas_daxpy(-ne*kappa, &E0.vector, &dE0.vector);
				if ( status != GSL_SUCCESS ) return status;
			}
			else {
				// select rows from the dE and E matrices
				gsl_vector_view dEi = gsl_matrix_row(&dE_vecs.matrix, i);
				gsl_vector_const_view Enext = gsl_matrix_const_row(&E_vecs.matrix, i);
				gsl_vector_const_view Eprev = gsl_matrix_const_row(&E_vecs.matrix, i-1);
				// compute dE
				status = gsl_blas_daxpy(ne*kappa, &Eprev.vector, &dEi.vector); // y -> y + a x
				if ( status != GSL_SUCCESS ) return status;
				status = gsl_blas_daxpy(-ne*kappa, &Enext.vector, &dEi.vector);
				if ( status != GSL_SUCCESS ) return status;
			}
		}
		/* Compute the first dI vector (using the last E vector)
		 * dI_1/dt = ne*kappa E_ne - ni*gamma I_1
		 */
		gsl_vector_view dI0 = gsl_matrix_row(&dI_vecs.matrix, 0);
		gsl_vector_const_view I0 = gsl_matrix_const_row(&I_vecs.matrix, 0);
		gsl_vector_const_view Elast = gsl_matrix_const_row(&E_vecs.matrix, ne-1);
		// compute dI
		status = gsl_blas_daxpy(ne*kappa, &Elast.vector, &dI0.vector); // y -> y + a x
		if ( status != GSL_SUCCESS ) return status;
		status = gsl_blas_daxpy(-ni*gamma, &I0.vector, &dI0.vector);
		if ( status != GSL_SUCCESS ) return status;						
	} // if ( ne > 0 ), i.e. there is at least one E-compartment
	else { // else: n = 0
		/* we have a special case when ne = 0:
		 * dI_1/dt = -dS/dt - ni*gamma I_1
		 */
		gsl_vector_view dI0 = gsl_matrix_row(&dI_vecs.matrix, 0);
		gsl_vector_const_view I0 = gsl_matrix_const_row(&I_vecs.matrix, 0);
		// compute dI
		status = gsl_blas_daxpy(-1.0, &dS_vec.vector, &dI0.vector); // y -> y + a x
		if ( status != GSL_SUCCESS ) return status;
		status = gsl_blas_daxpy(-ni*gamma, &I0.vector, &dI0.vector);
		if ( status != GSL_SUCCESS ) return status;
	} // if/else ne >/= 0
	
	/* Above we computed:
	 * dI_1/dt = -dS/dt - ni*gamma I_1 (if ne = 0)
	 * dI_1/dt = ne*kappa E_ne - ni*gamma I_1 (if ne > 0)
	 * Now, we have to finish with
	 * dI_2/dt = ni*gamma I_1 - ni*gamma I_2
	 * ...
	 */
	for ( int i = 1; i < ni; ++i ) { // not that i starts at 1
		// select rows from the dI and I matrices
		gsl_vector_view dIi = gsl_matrix_row(&dI_vecs.matrix, i);
		gsl_vector_const_view Iprev = gsl_matrix_const_row(&I_vecs.matrix, i-1);
		gsl_vector_const_view Inext = gsl_matrix_const_row(&I_vecs.matrix, i);
		// compute dI
		status = gsl_blas_daxpy(ni*gamma, &Iprev.vector, &dIi.vector); // y -> y + a x
		if ( status != GSL_SUCCESS ) return status;
		status = gsl_blas_daxpy(-ni*gamma, &Inext.vector, &dIi.vector);
		if ( status != GSL_SUCCESS ) return status;
	}
			
	gsl_vector_free(I_sum);
	
	return status;
}

int submatrix_diagonal_add_constant(gsl_matrix* mat, int i, int j, int n, double x) {
	/* add x * I \in R^{n\times n} to the sub-matrix 
	 * starting at the (i,j)-th element of mat.
	 * 
	 *            [0 0 0]
	 * e.g. mat = [0 0 0]
	 *            [0 0 0]
	 * 
	 * x = 1, i = 1, j = 0, n = 2 results in 
	 * 
	 *       [0 0 0]
	 * mat = [1 0 0]
	 *       [0 1 0]
	 */
	gsl_matrix_view submat = gsl_matrix_submatrix(mat, i, j, n, n);
	gsl_vector_view diag_submat = gsl_matrix_diagonal(&submat.matrix);
	return gsl_vector_add_constant(&diag_submat.vector, x);

}

int initial_general_model(double epsilon, const OdeParamStruct & pars, double* y0) {
	/* don't use the full jacobian, but only the part for the infected
	 * classes to reduce computation time, but also to get rid of the
	 * egenvalue 0 that would occur for pertubations that only involve S.
	 * In this example, ne = 2 and ni = 2.
	 * 
	 *     [ 0  0   0  -B  -B ]
	 * J = [ 0 -K   0   B   B ]
	 *     [ 0  K  -K   0   0 ]
	 *     [ 0  0   K  -D   0 ]
	 *     [ 0  0   0   D  -D ]
	 * 
	 * with B = beta * diag(S) * C, D = ni * gamma * 1, K = ne * kappa * 1
	 * 
	 * j = [-K  0  B  B ] 
	 *     [ K -K  0  0 ]
	 *     [ 0  K -D  0 ]
	 *     [ 0  0  D -D ]
	 * 
	 * let s = -B/r (i1 + i2), then [s;e1;e2;i1;i2] is an eigenvector
	 * of J with eigenvalue r, if [e1;e2;i1;i2] is an eigenvector of
	 * j with eigenvalue r
	 */
	 
 	int status = GSL_SUCCESS; // return value TODO: use below...
	
	double* C = pars.C;
	double* S0 = pars.S0;
	double beta = pars.beta;
	double gamma = pars.gamma;
	double kappa = pars.kappa;
	int ne = pars.ne;
	int ni = pars.ni;
	int na = pars.na;
	
	gsl_matrix_const_view C_mat = gsl_matrix_const_view_array(C, na, na);
	gsl_vector_const_view S0_vec = gsl_vector_const_view_array(S0, na);
	gsl_matrix* diag_S0_mat = gsl_matrix_calloc(na, na); // ALLOC
	gsl_vector_view diag_S0_vec = gsl_matrix_diagonal(diag_S0_mat);
	gsl_vector_memcpy(&diag_S0_vec.vector, &S0_vec.vector);
		
	// construct jacobian
	gsl_matrix* jac = gsl_matrix_calloc((ne+ni)*na, (ne+ni)*na); // calloc sets values to 0.0
	
	for ( int i = 0; i < ni; ++i ) { // top row, right blocks
		gsl_matrix_view subjac = gsl_matrix_submatrix(jac, 0*na, (ne+i)*na, na, na);
		if ( i == 0 ) {
			// compute B the first time
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, beta, diag_S0_mat, &C_mat.matrix, 1.0, &subjac.matrix);
		}
		else {
			// copy the first B for the others
			gsl_matrix_const_view subjac0 = gsl_matrix_const_submatrix(jac, 0*na, ne*na, na, na);
			gsl_matrix_add(&subjac.matrix, &subjac0.matrix);
		}
		// notice that the jac is initiated zero
	}
	for ( int i = 0; i < ne; ++i ) { // E-part of the diagonal
		submatrix_diagonal_add_constant(jac, i*na, i*na, na, -ne*kappa);
	}
	for ( int i = 0; i < ni; ++i ) { // I-part of the diagonal
		submatrix_diagonal_add_constant(jac, (ne+i)*na, (ne+i)*na, na, -ni*gamma);
	}
	for ( int i = 0; i < ne; ++i ) { // E-part of the sub-diagonal
		submatrix_diagonal_add_constant(jac, (i+1)*na, i*na, na, ne*kappa);
	}
	for ( int i = 0; i < ni-1; ++i ) { // I-part of the sub-diagonal
		submatrix_diagonal_add_constant(jac, (i+1+ne)*na, (i+ne)*na, na, ni*gamma);
	}
	
	// compute eigenspace and select right eigenvector
	gsl_vector_complex* eval = gsl_vector_complex_alloc((ne+ni)*na); // ALLOC
	gsl_matrix_complex* evec = gsl_matrix_complex_alloc((ne+ni)*na, (ne+ni)*na); // MALLOC
	gsl_eigen_nonsymmv_workspace* w = gsl_eigen_nonsymmv_alloc((ne+ni)*na); // MALLOC
	gsl_eigen_nonsymmv(jac, eval, evec, w);
	gsl_vector_view re_eval = gsl_vector_complex_real(eval);
	size_t idx = gsl_vector_max_index(&re_eval.vector);
	double r0 = gsl_vector_get(&re_eval.vector, idx);
	gsl_vector_complex_view h_complex = gsl_matrix_complex_column(evec, idx);
	gsl_vector_view h = gsl_vector_complex_real(&h_complex.vector);
	
	// store the initial state in y0
	gsl_vector* I_sum = gsl_vector_calloc(na); // ALLOC
	gsl_vector_view y0_vec = gsl_vector_view_array(y0, (1+ne+ni)*na);
	gsl_vector_set_zero(&y0_vec.vector);
	gsl_vector_view IE_vec = gsl_vector_subvector(&y0_vec.vector, na, (ne+ni)*na);
	gsl_vector_memcpy(&IE_vec.vector, &h.vector);
	// and create s by s = -B * (i1+i2+...) / r0
	gsl_vector_view S_vec = gsl_vector_subvector(&y0_vec.vector, 0, na);
	for ( int i = 0; i < ni; ++i ) {
		gsl_vector_const_view Ii_vec = gsl_vector_const_subvector(&IE_vec.vector, (ne+i)*na, na);
		gsl_vector_add(I_sum, &Ii_vec.vector); 
	}
	if ( r0 != 0.0 ) {
		// B was already computed: top right block of jac...
		gsl_matrix_const_view subjac = gsl_matrix_const_submatrix(jac, 0*na, (ne+ni-1)*na, na, na);
		gsl_blas_dgemv(CblasNoTrans, -1.0/r0, &subjac.matrix, I_sum, 0.0, &S_vec.vector);
	}
	else {
		std::cerr << "# WARNING: prevented division by zero while making"
		          << "initial condition" << RIGHT_HERE << std::endl;
	}

	// re-scale and test
	double norm = gsl_blas_dnrm2(&y0_vec.vector);
	gsl_vector_scale(&y0_vec.vector, 1.0/norm); // normalize
	if ( !gsl_vector_isnonneg(&IE_vec.vector) ) gsl_vector_scale(&y0_vec.vector, -1.0); // 'flip' the vector
	if ( !gsl_vector_isnonneg(&IE_vec.vector) ) { 
        // flipping was not sufficiant: mixed pos and neg terms (should not happen)...
		std::cerr << "# WARNING: cannot find a proper initial condition using the jacobian. "
		          << "using uniform vector in stead" << RIGHT_HERE << std::endl;
		// fall back on uniform initial condition
		gsl_vector_set_zero(&y0_vec.vector);
		gsl_vector_set_all(&S_vec.vector, -1.0);
		gsl_vector_set_all(&IE_vec.vector, 1.0/(ne+ni));
	}
	gsl_vector_scale(&y0_vec.vector, epsilon); // make the deviation small
	gsl_vector_add(&S_vec.vector, &S0_vec.vector); // add the disease-free steady state
	if ( !gsl_vector_isnonneg(&y0_vec.vector) ) { // should not happen
		std::cerr << "# WARNING: initial condition with negative values. "
		          << "replacing negative values with 0.0" << RIGHT_HERE << std::endl;
		// remove negative elements
		for ( int i=0; i < (1+ne+ni)*na; ++i ) y0[i] = ( y0[i] >= 0.0 ? y0[i] : 0.0 );
	}
	// free jacobian matrix and aux matrices, vectors and memory
	gsl_matrix_free(jac); // FREE
	gsl_matrix_free(diag_S0_mat); // FREE
	gsl_eigen_nonsymmv_free(w); // FREE
	gsl_vector_complex_free(eval); // FREE
	gsl_matrix_complex_free(evec); // FREE
	gsl_vector_free(I_sum); // FREE
	
	return status;
}

bool integrate_general_model(const double ts[], double** ss, int len, OdeParamStruct & pars) {
	double epsilon = INOCULUM_SIZE; // used for scaling pertubation from disease free equilibrium.
	int na = pars.na; // number of age classes
	int ne = pars.ne; // number of exposes stages
	int ni = pars.ni; // number of infectious stages
	double t0 = pars.t0;
	double* S0 = pars.S0;
	
	gsl_vector_const_view S0_vec = gsl_vector_const_view_array(S0, na);

	size_t dim = (1+ne+ni)*na;
	gsl_odeiv2_system sys = {dy_general_model, NULL, dim, &pars};
	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0); // ALLOC
	int status = GSL_SUCCESS;

	double t = t0;
	double* y = new double[(1+ne+ni)*na]; // NEW
	gsl_vector_view S_vec = gsl_vector_view_array(y, na); // for copying to ss
	status = initial_general_model(epsilon, pars, y); // initiate y
	
	if ( status != GSL_SUCCESS ) {
		std::cout << "# WARNING! bad initial condition status" << RIGHT_HERE << std::endl;
	}
	
	gsl_vector_view ss_vec; // to be used below...
	for ( int i = 0; i < len; ++i ) {
		if ( t0 >= ts[i] ) { // pre-epidemic
			ss_vec = gsl_vector_view_array(ss[i], na);
			gsl_vector_memcpy(&ss_vec.vector, &S0_vec.vector);
		}
		else { // during epidemic (todo: stop when prevalence below threshold)
			status = gsl_odeiv2_driver_apply(d, &t, ts[i], y); // y gets updated...
			
			if ( status != GSL_SUCCESS ) {
				std::cout << "# WARNING! bad integrator status" << RIGHT_HERE << std::endl;
				break;
			}
			ss_vec = gsl_vector_view_array(ss[i], na);
			gsl_vector_memcpy(&ss_vec.vector, &S_vec.vector);
		}
	}
	gsl_odeiv2_driver_free(d); // FREE
	delete[] y; // DELETE
	return (status == GSL_SUCCESS);
}


double calcR0(const OdeParamStruct & pars) {
	// todo: remove calcR0 from aux.cpp
	return calcR0(pars.S0, pars.C, pars.na, pars.beta, pars.gamma);
}
