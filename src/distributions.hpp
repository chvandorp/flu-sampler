#ifndef DISTRIBUTIONS_HPP
#define DISTRIBUTIONS_HPP

#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>


/* simply use the poisson distribution for dispersion (d) smaller than 
 * the following lower bound
 */
#define POLYA_DISPERSION_LOWER_BOUND 1e-4
#define LOG2PI 1.8378770664093 // redundant: same as M_LNPI + M_LN2

// a very small wrapper for gls_rng
class Rng {
public:
	Rng(); // alloc the gsl_rng, std seed
	Rng(unsigned long ); // alloc the gsl_rng, seed
	~Rng(); // free the rng
	void seed(unsigned long );
	// discete distributions
	unsigned long integer();
	unsigned long integer(unsigned long ); // mod argument
	bool bernoulli(double ); // p
	int binomial(int , double ); // n, p
	int poisson(double ); // mean
	int negBinomial(double , double ); // mean, n (1/dispersion)
	int polya(double , double ); // mean, d (dispersion)
	// continuous disributions
	double uniform(); // Unifirm(0,1), excludes boundary
	double stdNormal(); // N(0,1)
	double bactrian(double ); // m, mixture of two Gaussians. 0 <= m < 1
	gsl_rng* r; // todo: make private, adjust carnes.hpp/cpp
private:
};

// likelihoods

double ll_binomial(int , double , int ); // n, p, k
double ll_poisson(double , int ); // lambda, k
double ll_neg_binomial(double , double , int ); // lambda=np/(1-p), n, k
double ll_polya(double , double , int ); // lambda=np/(1-p), d=1/n, k

// Priors

class Prior {
public:
	Prior();
	virtual ~Prior();
	virtual Prior* dup() const; // make a copy (borrowed from the D language)
	virtual double loglike(double ) const; // returns 0.0
protected:	
};

class UniformPrior : public Prior {
public:
	UniformPrior(double , double );
	Prior* dup() const; // make a copy
	double loglike(double ) const;
protected:
	double a, b;
	double ll;
};

class BetaPrior : public Prior {
public:
	BetaPrior(double , double );
	Prior* dup() const; // make a copy
	double loglike(double ) const;
protected:
	double alpha, beta, logBeta;
};

class NormalPrior : public Prior {
public:
	NormalPrior(double , double ); // mu, sigmasq
	Prior* dup() const; // make a copy
	double loglike(double ) const;
protected:
	double mu, sigmasq, logsigmasq;
};

class StdNormalPrior : public Prior {
public:
	StdNormalPrior();
	Prior* dup() const;
	double loglike(double ) const;
protected:
};

class HalfNormalPrior : public Prior {
public:
	HalfNormalPrior(double ); // sigmasq
	Prior* dup() const; // make a copy
	double loglike(double ) const;
protected:
	double sigmasq; // variance of the 'full' normal distribution
	double logsigmasq; 
};

class LogStdNormalPrior : public Prior {
	LogStdNormalPrior();
	Prior* dup() const; // make a copy
	double loglike(double ) const;
protected:
};

class GammaPrior : public Prior {
public:
	GammaPrior(double , double ); // rate, shape
	Prior* dup() const;
	double loglike(double ) const;
protected:
	double rate, shape; // expectation is shape/rate
	double lograte, loggammashape; // pre-compute once!!
};

class InvGammaPrior : public Prior {
public:
	InvGammaPrior(double , double ); // rate, shape
	Prior* dup() const; // make a copy
	double loglike(double ) const;
protected:
	double rate, shape;
	double lograte, loggammashape; // pre-compute once!!
};

// Transformations Used for hypo-priors

class Transformation { // the base class acts as the identity
public:
	Transformation();
	virtual ~Transformation();
	virtual double evalFun(double );
	virtual double evalJac(double );
	virtual double evalLogJac(double );
protected:	
};

class LinearTransformation : public Transformation {
public:
	LinearTransformation(double , double );
	double evalFun(double );
	double evalJac(double );
	double evalLogJac(double );
protected:
	double loc, scale;
};

class LogitTransformation : public Transformation {
public:
	LogitTransformation(double , double );
	double evalFun(double );
	double evalJac(double );
	double evalLogJac(double );
protected:
	double loc, scale;
};

class LogTransformation : public Transformation {
public:
	LogTransformation(double , double );
	double evalFun(double );
	double evalJac(double );
	double evalLogJac(double );
protected:
	double loc, scale;	
};

#endif
