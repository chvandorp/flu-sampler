#ifndef MCMC_STATE_HPP
#define MCMC_STATE_HPP

#include <iostream>
#include <cmath>
#include <list>
#include <map>
#include <algorithm> // std::advance
#include <regex> // for finding parameter names

// todo: these files should be included by mcmc_full_state.hpp and/or mcmc_season_state.hpp

#include "data.hpp"
#include "integrator.hpp"
#include "distributions.hpp"
#include "parameters.hpp"
#include "locks.hpp"
#include "constants.hpp"
#include "aux.hpp"
#include "parallelism.hpp"
#include "exceptions.hpp"

enum McmcMode {	WARMUP_MODE=0, SAMPLING_MODE };

// todo: further implement a parent class for McmcFullState and McmcSeasonState
class McmcState {
public:
	McmcState(); // default constructor
	McmcState(const McmcState & ); // copy constructor
	McmcState & operator=(const McmcState & ); // copy assignment operator
	virtual ~McmcState(); // calls the private clear method
	virtual McmcState* dup() const; // returns a copy of self, override in derived classes!
	// likelihood methods
	virtual double loglike(); // log P(D|theta)
	virtual double prior_loglike() const; // log pi(theta)
	/* log pi(theta) * P(D|theta), equals P(theta|D) up to a constant 
	 * the argument is the 'temperature' used for path sampling.
	 * non-const, because it sets the ll_cache
	 */
	virtual double posterior_loglike(double );
	// general methods (todo: make set cache functions protected?)
	bool getLoglikeCache(double & ) const; // returns ll_cache_good
	void setLoglikeCache(double ); // assumes that the argument is correct...
	void resetLoglikeCache(); // goes back to previous ll

	bool getPostLoglikeCache(double & ) const; // returns pll_cache_good
	void setPostLoglikeCache(double ); // assumes that the argument is correct...
	void resetPostLoglikeCache(); // goes back to previous pll

	void setLocks(const std::list<Lock> & ); // use a list of "locks" to restrict the model
	virtual int countParameters() const; // count all unlocked parameters
	void setMode(McmcMode ); // set mode, also for submodels
	std::string getName() const;
	bool findValue(const std::string & , double & ) const; 
	// find the value of the parameter with a certain name
protected:
	std::string name; // identifier
	// parameters
	std::map<std::string, ParPrior*> parameters;
	// submodels
	std::map<std::string, McmcState*> submodels;
	// caches
	double ll_cache; // cache for the likelihood function log p(D|theta)
	double old_ll; // used for by resetLoglikeCache()
	bool ll_cache_good;
	// save pll for MAP estimates (TODO)
	double pll_cache;
	double old_pll; // FIXME: these is also a variable in use called pll_old (in updateSinglePar)...
	bool pll_cache_good;
	McmcMode mode;
private:
	virtual void copy(const McmcState & );
	virtual void clear();
	virtual void mapParameters(); // should at some point be removed
};

#endif
