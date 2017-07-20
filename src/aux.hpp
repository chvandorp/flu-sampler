#ifndef AUX_HPP
#define AUX_HPP

#include <iostream> // testing
#include <cmath> // log, exp
#include <algorithm> // max_element
#include <list>
#include <string>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include "exceptions.hpp"
#include "constants.hpp" // DAYS_IN_YEAR, M_PI, YEAR_ZERO

double logit(double );
double invLogit(double );
double calcTrendendProb(double , double , int , double c=0.0); // centralize around c
double calcLinearTrend(double , double , int , double c=0.0); // centralize around c
double calcBackgroundIli(double , double , double , int );
double calcR0(double* , double* , int , double , double ); // S0, C, size, beta, gamma

double calc_elpd(const std::list<double> & );
double mean(const std::list<double> & );
double sample_variance(const std::list<double> & );
double log_sum_exp(const std::list<double> & );
double log_mean_exp(const std::list<double> & ); // log_sum_exp(xs) - log(xs.size())

double softpos(double x); // log(exp(x) + 1) (overflow protected)

int yearToSeason(int s); // Map the year y to the index of season y/y+1, e.g. 1969 -> 0 (69/70)

std::string seasonString(int s); // 0 -> 69/70

#endif
