#include "distributions.hpp"

// a wrapper for gsl_rng

Rng::Rng() {
	r = gsl_rng_alloc(gsl_rng_mt19937);
	seed(144169);
}
Rng::Rng(unsigned long s) {
	r = gsl_rng_alloc(gsl_rng_mt19937);
	seed(s);
}
Rng::~Rng() { 
	gsl_rng_free(r); 
}
void Rng::seed(unsigned long s) {
	gsl_rng_set(r, s);
}
unsigned long Rng::integer() {
	return gsl_rng_get(r);
}
unsigned long Rng::integer(unsigned long n) {
	return gsl_rng_uniform_int(r, n);
}
bool Rng::bernoulli(double p) {
	return gsl_ran_bernoulli(r, p);
}
int Rng::binomial(int n, double p) {
	return gsl_ran_binomial(r, p, n);
}
int Rng::poisson(double lambda) {
	return gsl_ran_poisson(r, lambda);
}
int Rng::negBinomial(double lambda, double n) {
	/* lambda is the mean = n p/(1-p), so
	 * p = lambda / (n + lambda)
	 */
	double p = lambda / (n+lambda);
	return gsl_ran_negative_binomial(r, 1.0-p, n); 
	// notice the 1-p, different from standard definition 
}
int Rng::polya(double lambda, double d) {
	if ( d < POLYA_DISPERSION_LOWER_BOUND ) {
		return poisson(lambda);
	}
	else {
		return negBinomial(lambda, 1/d);
	}
}
double Rng::uniform() {
	return gsl_rng_uniform_pos(r);
}
double Rng::stdNormal() {
	return gsl_ran_gaussian(r, 1.0);
}
double Rng::bactrian(double m) { 
	/* for proposals in the shape of a Bactrian Camel /\/\
	 * see Yang and Rodriguez, PNAS (2013). when m=0, this
	 * is simply the standard normal distribution
	 */
	 return sqrt(1-m*m) * gsl_ran_gaussian(r, 1.0) + ( bernoulli(0.5) ? m : -m );
}


// likelihoods...

double ll_binomial(int n, double p, int k) {
	/* likelihood = {n choose k} * p^k * (1-p)^{n-k}
	 */
	return  gsl_sf_lnchoose(n, k) + k * log(p) + (n-k) * log(1-p);
}

double ll_poisson(double lambda, int k) {
	/* likelihood = lambda^k/k! * exp(-lambda), 
	 * so log likelihood = k * log(lambda) - log(k!) - lambda
	 */
	return k * log(lambda) - lambda - gsl_sf_lngamma(k+1);
}

double ll_neg_binomial(double lambda, double n, int k) {
	/* likelihood is Gamma(n + k) / (Gamma(k+1) Gamma(n)) * p^n * (1-p)^k
	 * where p = lambda / (n + lambda) and n = 1/d
	 */
	double p = lambda / (n+lambda);
	return gsl_sf_lngamma(n+k) - gsl_sf_lngamma(k+1) - gsl_sf_lngamma(n) + 
		n * log(1.0-p) + k * log(p);
}

double ll_polya(double lambda, double d, int k) {
	if ( d < POLYA_DISPERSION_LOWER_BOUND ) {
		return ll_poisson(lambda, k);
	}
	else {
		return ll_neg_binomial(lambda, 1/d, k);
	}
}


/* methods for priors */

Prior::Prior() { /* empty */ }
Prior::~Prior() { /* empty */ }
Prior* Prior::dup() const {
	return new Prior(*this);
}
double Prior::loglike(double x) const {
	return 0.0;
}

UniformPrior::UniformPrior(double a, double b) : a(a), b(b) {
	ll = -log(b - a); // todo: throw something when ill defined
}
Prior* UniformPrior::dup() const {
	return new UniformPrior(*this);
}
double UniformPrior::loglike(double x) const {
	return ll; // NB: does not check that x is between a and b... todo
}

BetaPrior::BetaPrior(double alpha, double beta) : alpha(alpha), beta(beta) {
	logBeta = gsl_sf_lnbeta(alpha, beta);
}
Prior* BetaPrior::dup() const {
	return new BetaPrior(*this);
}
double BetaPrior::loglike(double x) const {
	return (alpha-1.0)*log(x) + (beta-1.0)*log(1.0-x) - logBeta;
}

NormalPrior::NormalPrior(double mu, double sigmasq) : mu(mu), sigmasq(sigmasq) {
	logsigmasq = log(sigmasq);
}
Prior* NormalPrior::dup() const {
	return new NormalPrior(*this);
}
double NormalPrior::loglike(double x) const {
	// return -0.5*(log(2*M_PI) + logsigmasq) - (x-mu)*(x-mu)/(2.0*sigmasq);
	return -0.5*(M_LN2 + M_LNPI + logsigmasq) - (x-mu)*(x-mu)/(2.0*sigmasq);
}

HalfNormalPrior::HalfNormalPrior(double sigmasq) : sigmasq(sigmasq) {
	logsigmasq = log(sigmasq);
}
Prior* HalfNormalPrior::dup() const {
	return new HalfNormalPrior(*this);
}
double HalfNormalPrior::loglike(double x) const {
	// return -0.5*log(0.5*M_PI*sigmasq) - x*x/(2.0*sigmasq);
	return -0.5*(M_LNPI - M_LN2 + logsigmasq) - x*x/(2.0*sigmasq);
}

StdNormalPrior::StdNormalPrior() {
	/* empty */
}
Prior* StdNormalPrior::dup() const {
	return new StdNormalPrior();
}
double StdNormalPrior::loglike(double x) const {
	// return -0.5*LOG2PI - x*x/2.0;
	return -0.5*(M_LNPI + M_LN2) - x*x/2.0;
}

LogStdNormalPrior::LogStdNormalPrior() {
	/* empty */
}
Prior* LogStdNormalPrior::dup() const {
	return new LogStdNormalPrior();
}
double LogStdNormalPrior::loglike(double x) const {
	double logx = log(x);
	return -0.5*(M_LNPI + M_LN2) - logx - logx*logx/2.0; // like is 1/(x * sqrt(2pi)) * exp(-log(x)^2/2)
}

GammaPrior::GammaPrior(double rate, double shape) : rate(rate), shape(shape) {
	lograte = log(rate);
	loggammashape = gsl_sf_lngamma(shape);	
}
Prior* GammaPrior::dup() const {
	return new GammaPrior(*this);
}
double GammaPrior::loglike(double x) const {
	return shape*lograte - loggammashape + (shape-1)*log(x) - rate*x;
}


InvGammaPrior::InvGammaPrior(double rate, double shape) : rate(rate), shape(shape) {
	lograte = log(rate);
	loggammashape = gsl_sf_lngamma(shape);
}
Prior* InvGammaPrior::dup() const {
	return new InvGammaPrior(*this);
}
double InvGammaPrior::loglike(double x) const {
	return shape*lograte - loggammashape + (-shape-1)*log(x) - rate/x;
}

/* methods for transformations */

Transformation::Transformation() {
	/* empty */
}
Transformation::~Transformation() {
	/* empty */
}
double Transformation::evalFun(double x) {
	return x;
}
double Transformation::evalJac(double x) {
	return 1.0;
}
double Transformation::evalLogJac(double x) {
	return 0.0;
}


LinearTransformation::LinearTransformation(double loc, double scale) :
		loc(loc), scale(scale) {
	/* empty */		
}
double LinearTransformation::evalFun(double x) {
	return (x - loc) / scale;
}
double LinearTransformation::evalJac(double x) {
	return 1.0/scale;
}
double LinearTransformation::evalLogJac(double x) {
	return -log(scale);
}


LogitTransformation::LogitTransformation(double loc, double scale) :
		loc(loc), scale(scale) {
	/* empty */		
}
double LogitTransformation::evalFun(double x) {
	return (log(x/(1.0-x)) - loc) / scale;
}
double LogitTransformation::evalJac(double x) {
	return 1.0/(scale*x*(1.0-x));
}
double LogitTransformation::evalLogJac(double x) {
	return -log(scale*x*(1.0-x));
}


LogTransformation::LogTransformation(double loc, double scale) :
		loc(loc), scale(scale) {
	/* empty */
}
double LogTransformation::evalFun(double x) {
	return (log(x) - loc) / scale;
}
double LogTransformation::evalJac(double x) {
	return 1/(x*scale);
}
double LogTransformation::evalLogJac(double x) {
	return -log(x*scale);
}
