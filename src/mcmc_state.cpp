#include "mcmc_state.hpp"

McmcState::McmcState() : ll_cache(0.0), old_ll(0.0), ll_cache_good(false), 
			pll_cache(0.0), old_pll(0.0), pll_cache_good(false), 
			mode(WARMUP_MODE) { /* empty */ }
McmcState::McmcState(const McmcState & state) {
	copy(state);
}
McmcState & McmcState::operator=(const McmcState & state) {
	if ( this != &state ) {
		clear();
		copy(state);
	}
	return *this;
}
McmcState::~McmcState() {
	clear();
}
McmcState* McmcState::dup() const {
	return new McmcState(*this);
}
void McmcState::copy(const McmcState & state) {
	name = state.name;
	ll_cache = state.ll_cache;
	old_ll = state.old_ll;
	ll_cache_good = state.ll_cache_good;
	pll_cache = state.pll_cache;
	old_pll = state.old_pll;
	pll_cache_good = state.pll_cache_good;
	mode = state.mode;
	
	// in the future, copy the parameter and submodel maps
	/** TODO: copy the parameter map...
	for ( auto it = state.parameters.begin(); it != state.parameters.end(); ++it ) {
		parameters[it->first] = it->second->dup(); // dup makes a copy...
	}
	*/
	/** TODO: copy the submodels...
	for ( auto it = state.submodels.begin(); it != state.submodels.end(); ++it ) {
		submodels[it->first] = it->second->dup(); // dup is virtual and derived classes override it
	}
	*/
}
void McmcState::clear() {
	// in the future: delete submodels and parameters in the maps
	/**
	for ( auto it = parameters.begin(); it != parameters.end(); ++it ) {
		delete it->second;
	}
	for ( auto it = submodels.begin(); it != submodels.end(); ++it ) {
		delete it->second;
	}
	*/
}

double McmcState::prior_loglike() const { // log pi(theta)
	double pi = 0.0;
	std::map<std::string, ParPrior*>::const_iterator pit = parameters.begin();
	for ( ; pit != parameters.end(); ++pit ) {
		pi += pit->second->loglike();
	}
	std::map<std::string, McmcState*>::const_iterator mit = submodels.begin();
	for ( ; mit != submodels.end(); ++mit ) {
		pi += mit->second->prior_loglike();
	}
	return pi;
}
double McmcState::loglike() {
	double ll = 0.0;
	setLoglikeCache(ll);
	return ll;
}
double McmcState::posterior_loglike(double T) { 
	// log pi(theta) * P(D|theta), equals P(theta|D) up to a constant
	double pll = 0;
	if ( T != 0 ) pll = prior_loglike() + T * loglike();
	else pll = prior_loglike();
	setPostLoglikeCache(pll);
	return pll;
}

bool McmcState::getLoglikeCache(double & ll) const {
	ll = ll_cache;
	return ll_cache_good;
}
void McmcState::setLoglikeCache(double ll) {
	old_ll = ll_cache;
	ll_cache = ll;
	ll_cache_good = true;
}
void McmcState::resetLoglikeCache() {
	ll_cache = old_ll;
	ll_cache_good = true;
}	

bool McmcState::getPostLoglikeCache(double & pll) const {
	pll = pll_cache;
	return pll_cache_good;
}
void McmcState::setPostLoglikeCache(double pll) {
	old_pll = pll_cache;
	pll_cache = pll;
	pll_cache_good = true;
}
void McmcState::resetPostLoglikeCache() {
	pll_cache = old_pll;
	pll_cache_good = true;
}	

void McmcState::setLocks(const std::list<Lock> & locks) {
	for ( auto lit = locks.begin(); lit != locks.end(); ++lit ) {
		std::regex re_name(lit->name); // name might contain wildcards
		for ( auto pit = parameters.begin(); pit != parameters.end(); ++pit ) {
			if ( std::regex_match(pit->first, re_name) ) {
				if ( pit->second != NULL ) {
					ParPrior* pp = pit->second;
					switch ( lit->type ) {
						case Lock::INVALID: {
							// todo: warning?
							break;
						}
						case Lock::CONSTANT: {
							// simply lock the parameter to the target value
							pp->lock(lit->target_value);
							break;
						}
						case Lock::EQUALS: {
							// find the other paremeter...
							std::regex re_target(lit->target_name);
							// target_name might still contain wildcards
							for ( auto tpit = parameters.begin(); tpit != parameters.end(); ++tpit ) {
								if ( std::regex_match(tpit->first, re_target) ) {
									ParPrior* tpp = tpit->second;
									if ( tpp != NULL ) {
										// don't lock parameters to themselves
										if ( tpp != pp ) pp->lock(tpp);
									}
									else {
										std::stringstream ss;
										ss << "ref. to parameter '" << tpit->first << "' is NULL";
										throw MsgException(ss.str() + RIGHT_HERE);
									}
								}
							} // loop over parameters to find target
							break;
						}
						default: {
							// todo: warning?
							break;
						}
					}
				}
				else {
					std::stringstream ss;
					ss << "ref. to parameter '" << pit->first << "' is NULL";
					throw MsgException(ss.str() + RIGHT_HERE);
				}
			} // if regular expression matches parameter name
		} // for loop over parameters
	} // for loop over locks
	// set locks of sub-models
	for ( auto it = submodels.begin(); it != submodels.end(); ++it ) {
		McmcState* sm = it->second;
		if ( sm != NULL ) {
			sm->setLocks(locks);
		}
		else {
			std::stringstream ss;
			ss << "ref. to submodel '" << it->first << "' is NULL";
			throw MsgException(ss.str() + RIGHT_HERE);
		}
	}
}

int McmcState::countParameters() const {
	int k = 0;
	std::map<std::string, ParPrior*>::const_iterator pit = parameters.begin();
	for ( ; pit != parameters.end(); ++pit ) {
		ParPrior* par = pit->second;
		if ( par != NULL ) {
			if ( !par->isLocked() ) ++k;
		}
		else {
			std::stringstream ss;
			ss << "ref. to parameter '" << pit->first << "' is NULL";
			throw MsgException(ss.str() + RIGHT_HERE);
		}
	}
	std::map<std::string, McmcState*>::const_iterator mit = submodels.begin();
	for ( ; mit != submodels.end(); ++mit ) {
		McmcState* mod = mit->second;
		if ( mod != NULL ) {
			k += mod->countParameters();
		}
		else {
			std::stringstream ss;
			ss << "ref. to submodel '" << mit->first << "' is NULL";
			throw MsgException(ss.str() + RIGHT_HERE);
		}
	}
	return k;
}
void McmcState::setMode(McmcMode mode) {
	this->mode = mode;
	auto it = submodels.begin();
	for ( ; it != submodels.end(); ++it ) {
		McmcState* sm = it->second;
		if ( sm != NULL ) sm->setMode(mode);
	}
}
std::string McmcState::getName() const { return name; }

bool McmcState::findValue(const std::string & parname, double & value) const {
	auto pit = parameters.find(parname);
	if ( pit == parameters.end() ) return false;
	else if ( pit->second != NULL ) {
		value = pit->second->getValue();
		// std::cout << parname << " = " << value << std::endl; // TESTING
		return true;
	}
	else {
		std::stringstream ss;
		ss << "ref. to parameter '" << pit->first << "' is NULL";
		throw MsgException(ss.str() + RIGHT_HERE);
	}
}

void McmcState::mapParameters() {
	std::cout << "this message should never be printed" << RIGHT_HERE << std::endl;
}
