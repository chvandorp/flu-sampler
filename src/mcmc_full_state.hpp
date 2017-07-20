#ifndef MCMC_FULL_STATE_HPP
#define MCMC_FULL_STATE_HPP

#include "mcmc_state.hpp"
#include "mcmc_season_state.hpp"

class McmcFullState : public McmcState {
public:
	McmcFullState();
	// pass data to get number of time points for each season
	McmcFullState(const std::vector<IliDataStruct> & );
	McmcFullState(const McmcFullState & ); // copy constructor
	McmcFullState & operator=(const McmcFullState & ); // copy assignment operator
	~McmcFullState(); // destructor
	virtual McmcState* dup() const; // override
	// loglike functions
	double prior_loglike() const;
	// serial method
	double loglike(double , // beta
			std::vector<IliDataStruct> & , // data
			bool reintegrate=true,
			bool record_lls=false);
	// parallel method
	double loglike(double , // beta
			std::vector<IliDataStruct> & , // data
			WorkerPool & , // worker pool
			bool reintegrate=true,
			bool record_lls=false);
	// serial method
	double posterior_loglike(double , // 'temperature' for path sampling;
			double , // beta
			std::vector<IliDataStruct> & , // data
			bool reintegrate=true);
	// parallel method
	double posterior_loglike(double , // 'temperature' for path sampling;
			double , // beta
			std::vector<IliDataStruct> & , // data
			WorkerPool & , // worker pool
			bool reintegrate=true);
	std::vector<IliDataStruct> simulate(double , // beta
			const std::vector<IliDataStruct> & , // data
			Rng & ); // random number generator
	// serial update
	void update(double , // T
			double , // beta,
			std::vector<IliDataStruct> & , // data
			Rng & ); // random number generator
	// parallel update
	void update(double , // temperature for path sampling
			double , // beta,
			std::vector<IliDataStruct> & , // data
			Rng & , WorkerPool & ); // random number generator, worker pool
	void print(std::ostream & ) const; // pretty printing
	std::string xmlString() const;
	int countParameters() const;
protected:
	// auxilliary methods:
	void setDefaultValues(); // default values for the parameters / priors / bounds
	// serial
	void updateSinglePar(double ,	double , // T, beta, 
			std::vector<IliDataStruct> & , // data 
			ParPrior & , Rng & , bool reintegrate=true);
	// parallel
	void updateSinglePar(double ,	double , // T, beta, 
			std::vector<IliDataStruct> & , // data 
			ParPrior & , Rng & , WorkerPool & , 
			bool reintegrate=true);
	// members:
	// parameters...
	McmcSeasonState seasonStates[NUMBER_OF_SEASONS]; // automatically makes deep copies!
	ParPrior meanLogGamma; // mean recovery rate; hyperparameter
	ParPrior sdLogGamma; // sd of recovery rare; hyperparameter
	ParPrior meanLogKappa; // mean incubation rate; hyperparameter
	ParPrior sdLogKappa; // sd of incubation rare; hyperparameter
	ParPrior backgrounds[NUMBER_OF_AGE_CLASSES];
	ParPrior slopeBackgrounds[NUMBER_OF_AGE_CLASSES];
	ParPrior amplitudeBackground;
	ParPrior phaseBackground;
	ParPrior slopeReportingRates[NUMBER_OF_AGE_CLASSES];
	ParPrior offsetReportingRates[NUMBER_OF_AGE_CLASSES];
	ParPrior dispersion; // = 0 -> poisson, > 0 -> Polya (Negative Binomial)
	// hyper parameters
	ParPrior meanLogitS0[NUMBER_OF_AGE_CLASSES]; // offset of mean, really...
	ParPrior sdLogitS0[NUMBER_OF_AGE_CLASSES];
	ParPrior slopeLogitS0[NUMBER_OF_AGE_CLASSES];
	ParPrior meanT0;
	ParPrior sdT0;
	ParPrior sdSeasonalBgEffect;
private:
	void copy(const McmcFullState & );
	void clear();
	virtual void mapParameters();
};

std::ostream & operator<<(std::ostream & , const McmcFullState & );

#endif
