#include "parameters.hpp"

/* methods for ParPrior */

ParPrior::ParPrior() :
	name(""), value(0.0), old_value(0.0), pvar(1.0),
	lboundbool(false), lbound(0.0),
	uboundbool(false), ubound(0.0),
	homotopy_support(CONTRACTIBLE),
	prior(NULL), 
	stepcounter(0), acceptcounter(0), ac(0), uc(1),
	locked(false), lockedTo(NULL) {
	/* empty */
}

ParPrior::ParPrior(double value, double pvar, Prior* prior, std::string name) :
	name(name), value(value), old_value(value), pvar(pvar), prior(prior) {
	lboundbool = false;
	uboundbool = false;
	homotopy_support = CONTRACTIBLE; // todo: infer from prior??
	stepcounter = 0;
	acceptcounter = 0;
	ac = 0;
	uc = 1; // used in rescaling of pvar ($\exp(\pm1/\sqrt(ux))$)
	locked = false;
	lockedTo = NULL;
}
ParPrior::ParPrior(const ParPrior & pp) :
	name(pp.name), value(pp.getValue()), old_value(value), pvar(pp.pvar),
	lboundbool(pp.lboundbool), lbound(pp.lbound),
	uboundbool(pp.uboundbool), ubound(pp.ubound),
	homotopy_support(pp.homotopy_support),
	prior(pp.prior),
	stepcounter(pp.stepcounter), acceptcounter(pp.acceptcounter),
	ac(pp.ac), uc(pp.uc),
	locked(pp.locked), lockedTo(NULL) {
	if ( pp.prior != NULL )	prior = pp.prior->dup();
	/* NB: after copying, lockedTo is always NULL, 
	 * and value is the last value returned by pp.getValue()
	 */
}
ParPrior & ParPrior::operator=(const ParPrior & pp) {
	if ( this != &pp ) {
		name = pp.name;
		value = pp.getValue();
		old_value = pp.old_value;
		pvar = pp.pvar;
		lboundbool = pp.lboundbool;
		lbound = pp.lbound;
		uboundbool = pp.uboundbool;
		ubound = pp.ubound;        
		stepcounter = pp.stepcounter;
		acceptcounter = pp.acceptcounter;
		ac = pp.ac;
		uc = pp.uc;
		if ( pp.prior != NULL ) prior = pp.prior->dup();
		else prior = NULL;
		locked = pp.locked;
		lockedTo = NULL; // NB: after copying, lockedTo is always NULL
		homotopy_support = pp.homotopy_support;
	}
	return *this;
}
ParPrior::~ParPrior() {
	delete prior;
}
ParPrior* ParPrior::dup() const {
	return new ParPrior(*this);
}

void ParPrior::lock(double cval) { 
	value = cval;
	lockedTo = NULL;
	locked = true; 
}
void ParPrior::lock(ParPrior* pp) {
	locked = true;
	if ( pp != this ) lockedTo = pp;
	else std::cerr << "WARNING! trying to lock ParPrior to *this (ignored)." << std::endl;
}
bool ParPrior::isLocked() const { return locked; }
void ParPrior::setBounds(double lbound, double ubound) {
	this->lbound = lbound;
	lboundbool = true;
	this->ubound = ubound;
	uboundbool = true;
}
void ParPrior::setLBound(double lbound) {
	this->lbound = lbound;
	lboundbool = true;
}
void ParPrior::setUBound(double ubound) {
	this->ubound = ubound;
	uboundbool = true;
}
bool ParPrior::isBounded() const {
    return lboundbool && uboundbool;
}
double ParPrior::getLengthInterval() const {
    if ( isBounded() ) {
        return ubound - lbound;
    }
    else {
        return std::numeric_limits<double>::infinity();
    }
}
void ParPrior::setHomotopy(HomotopyClass hc) {
	switch ( hc ) {
		case CONTRACTIBLE: {
			homotopy_support = CONTRACTIBLE;
			break;
		}
		case CIRCULAR: {
			if ( lboundbool && uboundbool ) {
				homotopy_support = CIRCULAR;
			}
			else {
				std::cerr << "# WARNING: can't make an unbounded support CIRCULAR. "
				          << "ParPrior::setHomotopy ignored." << std::endl;
			}
			break;
		}
		default: {
			std::cerr << "# WARNING: invalid HomotopyClass given." << std::endl;
			break;
		}
	}
}
double ParPrior::getValue() const {
	if ( locked && lockedTo != NULL ) return lockedTo->getValue();
	else return value;
}
double ParPrior::getAr() const {
	if ( stepcounter > 0 ) {
		return double(acceptcounter)/stepcounter;
	}
	else return 0.0;
}
double ParPrior::loglike(Transformation* fun) const { // by default fun = NULL
	if ( !locked ) {
		if ( prior != NULL ) {
			if ( fun == NULL ) {
				return prior->loglike(value);
			}
			else {
				return prior->loglike(fun->evalFun(value)) + fun->evalLogJac(value);
			}
		}
		else { // assume that prior is not informative
			if ( fun == NULL ) {
				return 0.0;
			}
			else {
				return fun->evalLogJac(value);
			}
		}
	}
	else { // the parameter is locked
		return 0.0; // constants don't contribute to the prior ll
	}
}

/* if lboundbool, mirror in lbound. if uboundbool, mirror in ubound
 * if both uboundbool and lboundbool, mirror in both bounds,
 * by using a combination of mirroring and modulus
 */
void ParPrior::mutate(Rng & rng) {
	if ( !locked ) {
		old_value = value;
		value = value + rng.stdNormal() * sqrt(pvar); // normal proposal
		// value = value + rng.bactrian(0.9) * sqrt(pvar); // bactrian proposal
		if ( lboundbool ) {
			if ( uboundbool ) { // the 'difficult' case
				double interval = (ubound-lbound); // length interval
				switch ( homotopy_support ) {
					case CONTRACTIBLE: { // apply double mirroring
						value = fmod(value-lbound, 2*interval); 
						// remainder after dividion by denominator
						if ( value < 0.0 ) value += 2*interval; 
						// fmod does not always return non-negative numbers
						if ( value > interval ) value = 2*interval - value; 
						// mirror
						value += lbound; 
						// translate back to the 'true' interval	
						break;
					}
					case CIRCULAR: { // value modulo interval
						value = fmod(value-lbound, interval);
						if ( value < 0.0 ) value += interval; 
						// fmod can return negative
						value += lbound; // translate back
						break;
					}
					default: {
						std::cerr << "# WARNING: invalid HomotopyClass." << std::endl;
						break;
					}
				}
			}
			else { // mirror in lbound?
				if ( value < lbound ) {
					value = 2*lbound - value;
				}
			}
		}
		else { // lboundbool == false!!
			if ( uboundbool ) { // mirror in ubound?
				if ( value > ubound ) {
					value = 2*ubound - value;
				}
			}
			// else: uboundbool == false && uboundbool == false, so value is valid!
		}
		stepcounter++;
	} // if ! locked... else do nothing
}
void ParPrior::ac1up() {
	if ( !locked ) {
		ac++; // used by updatePvar
		acceptcounter++; // used for diagnostics
	}
}
void ParPrior::go_back() { 
	// todo: maybe rename: reject; include acceptance counts 
	if ( !locked ) {
		value = old_value;
	}
}
void ParPrior::updatePvar() {
	if ( !locked ) {
		if ( stepcounter % PVAR_UPDATE_INTERVAL == 0 ) {
			double ar = double(ac) / PVAR_UPDATE_INTERVAL;
			// lower, or increase the proposal variance
			if ( ar < OPTIMAL_ACCEPTANCE_RATE ) {
				pvar *= exp(-1.0/sqrt(uc)); // jumps are too big!
			}
			else {
				pvar *= exp(1.0/sqrt(uc)); // jumps are too small!
			}
			// do some checks to prevent explosions
			if ( isBounded() ) {
                double lenival = getLengthInterval();
				pvar = std::min(pvar, lenival*lenival);
			}
            /* sometimes parameters get "stuck" during burnin.
             * when the pvar get too small, reset to a higher value...
             * FIXME: this is just a hack...
             */
            // if ( pvar < 1e-8 ) pvar = 0.01; // sqrt(pvar) < 1e-4
			// update/reset counters
			ac = 0; // reset acceptance counter
			uc++; // increase update counter
		}
	} // if ! locked... else do nothing
}

std::string ParPrior::getName() const { return name; }
std::string ParPrior::xmlString() const {
	std::stringstream ss;
	ss << "<param "
	   << "name='" << name << "' "
	   << "val='" << getValue() << "' " // use getValue, parameter may be locked
	   << "pvar='" << pvar << "' "
	   << "ar='" << getAr() << "' " // total acceptance ratio
	   << "lock='" << std::boolalpha << locked << "' />";	
	return ss.str();
}

void ParPrior::print(std::ostream & os) const {
	if ( !locked ) {
		// print pvar in a different color when acceptance rate is abnormal
		double ar = getAr();
		std::string font = DEF_TERM_FONT;
		if ( ar < LOW_ACCEPTANCE_RATE ) font = CYAN_TERM_FONT;
		else if ( ar > HIGH_ACCEPTANCE_RATE ) font = YELLOW_TERM_FONT;
		os << value << " (" << font << sqrt(pvar) << DEF_TERM_FONT << ")";
	}
	else {
		os << BBLUE_TERM_FONT << getValue() << DEF_TERM_FONT;
	}
}

std::ostream & operator<<(std::ostream & os, const ParPrior & pp) {
	pp.print(os);
	return os;
}


/** methods for CompundPar **/

CompoundPar & CompoundPar::operator=(double x) {
	setValue(x);
	return *this;
}
bool CompoundPar::operator<(double x) const { return getValue() < x; }
bool CompoundPar::operator>(double x) const { return getValue() > x; }
bool CompoundPar::operator==(double x) const { return getValue() == x; }
bool CompoundPar::operator<=(double x) const { 
	return *this < x || *this == x;
}
bool CompoundPar::operator>=(double x) const { 
	return *this > x || *this == x;
}
void CompoundPar::setValue(double value) { this->value = value; }
double CompoundPar::getValue() const { return value; }
std::string CompoundPar::xmlString() const {
	std::stringstream ss;
	ss << "<param name='" << name << "' "
	   << "val='" << value << "' />";
	return ss.str();
}
std::string CompoundPar::getName() const { return name; }
void CompoundPar::print(std::ostream & os) const {
	os << value;
}

std::ostream & operator<<(std::ostream & os, const CompoundPar & par) {
	par.print(os);
	return os;
}


/** Methods for DiscreteParPrior **/

DiscreteParPrior::DiscreteParPrior() : name(""), value(0), old_value(0), 
	jump_prob(1.0), lbound(0), ubound(0), 
	lboundbool(false), uboundbool(false) {
	/* empty */
}
DiscreteParPrior::DiscreteParPrior(int value, double jump_prob, std::string name) :
	name(name), value(value), old_value(value), jump_prob(jump_prob), 
	lbound(0), ubound(0), lboundbool(false), uboundbool(false) {
	/* empty */
}
void DiscreteParPrior::setBounds(int lbound, int ubound) {
	this->lbound = lbound;
	lboundbool = true;
	this->ubound = ubound;
	uboundbool = true;	
}
void DiscreteParPrior::setLBound(int lbound) { 
	this->lbound = lbound;
	lboundbool = true;
}
void DiscreteParPrior::setUBound(int ubound) {
	this->ubound = ubound;
	uboundbool = true;
}
int DiscreteParPrior::getValue() const {
	return value;
}
void DiscreteParPrior::mutate(Rng & rng) {
	if ( rng.bernoulli(jump_prob) ) {
		// head or tails
		if ( rng.bernoulli(0.5) || (lboundbool && value==lbound) ) ++value;
		else --value;
		// check bounds
	}
}
void DiscreteParPrior::go_back() {
	value = old_value;
}
void DiscreteParPrior::lock(int n) { 
	locked = true;
	value = n;
}
bool DiscreteParPrior::isLocked() const {
	return locked;
}

