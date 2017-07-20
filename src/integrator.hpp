#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include "aux.hpp" // calcR0 (todo: move calcR0 here completely)
#include "constants.hpp"
#include "exceptions.hpp"
 
struct OdeParamStruct {
	double* C; // contact matrix
	double* S0; // initial values susceptibles
	double t0; // start of the epidemic
	double beta; // infectionrate
	double gamma; // 1 / infectious
	double kappa; // 1 / exposed 
	int na; // number of age classes
	int ne; // number of exposed stages
	int ni; // number of infectious stages
};

bool integrate_general_model(const double[], double**, int, OdeParamStruct & );

double calcR0(const OdeParamStruct & );


#endif
