#ifndef MCMC_HPP
#define MCMC_HPP

#include <iostream>
#include <vector>
#include <list>

#include "data.hpp"
#include "constants.hpp"
#include "parallelism.hpp"
#include "eta.hpp" // for counting down remaining time
#include "mcmc_full_state.hpp" // the elements of the chain...

class Mcmc {
public:
	Mcmc(int , int , int , bool, unsigned long , std::string ); 
	// length, burnin, thinning, seed, version_identifier
	~Mcmc();
	void importData();
	void run();
	void runSingleSeason(int ); // for testing purposes
	void print(std::ostream & ) const;
	void simulateAndPrint(std::ostream & , int ); // can't be const: uses internal Rng
	double calcBic(); // stores result in bic_cache
	double calcWaic(); // stores result in waic_cache
	double calcWbic(); // stores result in wbic_cache
	void printPointwiseWaic(std::ostream & ) const;
protected:
	// auxilliary methods
	void saveMlState(const McmcFullState & ); // update both ml and mpl state
	// members
	int length;
	int burnin;
	int thinning;
	int print_thinning; // equal to thinning by default...
	bool wbic_mode; // sample at T=1/log(n) to calculate WBIC
	std::string vid; // version identifier (todo: rename to "name")
	std::list<McmcFullState> chain;
	McmcFullState mlState; // state with maximum likelihood
	McmcFullState mplState; // state with highest posterior likelihood
	// statistics
	double bic_cache; bool bic_cache_good;
	double waic_cache; bool waic_cache_good;
	double wbic_cache; bool wbic_cache_good;
	// data
	std::vector<IliDataStruct> data;
	double beta; // todo, make parameter (?)
	// tools
	Rng rng;
	WorkerPool workerPool;
};

std::ostream & operator<<(std::ostream & , const Mcmc & );

#endif
