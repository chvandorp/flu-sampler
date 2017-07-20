#include "aux.hpp"

double logit(double x) {
	return log(x/(1.0-x));
}
double invLogit(double x) {
	if ( x > 0.0 ) {
		return 1.0/(1.0+exp(-x));
	}
	else {
		double u = exp(x);
		return u/(u+1.0);
	}
}
double calcTrendendProb(double offset, double slope, int season, double c) {
	return invLogit(offset + slope * (season - c));
}
double calcLinearTrend(double offset, double slope, int season, double c) {
	return offset + slope * (season - c);
}
double calcBackgroundIli(double bg, double amp, double phase, int dayNr) {
	return invLogit(bg + amp*sin(2.0*M_PI*dayNr/DAYS_IN_YEAR + phase));
}

double calcR0(double* S0, double* C, int numAgeGrps, double beta, double gamma) {
	// make matrix 'views'
	gsl_matrix_const_view C_mat = gsl_matrix_const_view_array(C, numAgeGrps, numAgeGrps);
	gsl_vector_const_view S0_vec = gsl_vector_const_view_array(S0, numAgeGrps);

	// make 'diag(S0)'
    gsl_matrix* S0_mat = gsl_matrix_calloc(numAgeGrps, numAgeGrps);
    gsl_vector_view diag = gsl_matrix_diagonal(S0_mat);
    gsl_vector_memcpy(&diag.vector, &S0_vec.vector);
	
	// compute the product
	gsl_matrix* ngm = gsl_matrix_calloc(numAgeGrps, numAgeGrps);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, beta/gamma, S0_mat, &C_mat.matrix, 0.0, ngm);

	// compute the eigenvalues
	gsl_vector_complex* eval = gsl_vector_complex_alloc(numAgeGrps);
	gsl_matrix_complex* evec = gsl_matrix_complex_alloc(numAgeGrps, numAgeGrps);
	gsl_eigen_nonsymmv_workspace* w = gsl_eigen_nonsymmv_alloc(numAgeGrps);
	gsl_eigen_nonsymmv(ngm, eval, evec, w);
	gsl_eigen_nonsymmv_free(w);
	gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);
	
	// get the 'spectral radius'
	double R0 = gsl_complex_abs(gsl_vector_complex_get(eval, 0));
	
	// free stuff...
	gsl_matrix_free(ngm);
	gsl_matrix_free(S0_mat);
	gsl_vector_complex_free(eval);
	gsl_matrix_complex_free(evec);
		
	return R0;
}

double calc_elpd(const std::list<double> & lls) {
	if ( lls.size() < 2 ) {
		throw MsgException("can't compute elpd of loglike list of less than 2 elements" + RIGHT_HERE);
	}
	double lpd = log_mean_exp(lls);
	double pwaic = sample_variance(lls);
	return lpd - pwaic;
}

double log_sum_exp(const std::list<double> & xs) {
	// make a (non-const) copy that can be sorted...
	std::list<double> xs_prime = xs;
	xs_prime.sort();
	// find the maximum x to prevent overflow
	double xmax = *std::max_element(xs_prime.begin(), xs_prime.end());
	double z = 0.0;
	for ( auto it = xs_prime.begin(); it != xs_prime.end(); ++it ) {
		z += exp(*it - xmax);
	}
	return xmax + log(z);
}

double log_mean_exp(const std::list<double> & xs) {
	int n = xs.size();
	return log_sum_exp(xs) - log(n);
}

double mean(const std::list<double> & xs) {
	double z = 0.0;
	int n = xs.size();
	for ( auto it = xs.begin(); it != xs.end(); ++it ) {
		z += *it;
	}
	return z/n;
}

double sample_variance(const std::list<double> & xs) {
	// calc mean
	double m = mean(xs);
	double z = 0.0;
	int n = xs.size(); n--; // NB! sample variance
	for ( auto it = xs.begin(); it != xs.end(); ++it ) {
		z += (*it-m)*(*it-m);
	}
	return z/n;	
}

double softpos(double x) {
	if ( x > 0 ) return x + log(1 + exp(-x));
	else return log(exp(x) + 1);
}

int yearToSeason(int y) { 
	// Map the year y to the index of season y/y+1, e.g. 1969 -> 0 (69/70)
    return y - YEAR_ZERO + 1;
}


std::string seasonString(int s) {
    std::string str1 = std::to_string((YEAR_ZERO + s - 1) % 100);
    std::string str2 = std::to_string((YEAR_ZERO + s) % 100);
    if ( str1.size() < 2 ) str1 = "0" + str1; // padding of zeros
    if ( str2.size() < 2 ) str2 = "0" + str2;
    return str1 + "/" + str2;
}
