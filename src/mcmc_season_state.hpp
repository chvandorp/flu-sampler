#ifndef MCMC_SEASON_STATE_HPP
#define MCMC_SEASON_STATE_HPP

#include "mcmc_state.hpp"
#include <list>

class McmcSeasonState : public McmcState {
public:
	McmcSeasonState(); // trivial
	McmcSeasonState(int , int , int ); // argument is number of weeks, length of Ss, season number
	McmcSeasonState(const McmcSeasonState & );
	~McmcSeasonState(); // destructor: delete Ss, S0
	McmcSeasonState & operator=(const McmcSeasonState & );
	virtual McmcState* dup() const; // override
	double prior_loglike(const double[] , const double[] , // meanLogitS0, sdLogitS0
			double , double , // meanT0, sdT0
			double , double , // meanLogGamma, sdLogGamma
			double , double , // meanLogKappa, sdLogKappa
			double ) const; // sdSeasonalBgEffect
	double loglike(const double[] , const double[] , double , // background, reportingrate, dispersion
			double , // beta
			double , double , // oscilation of bg ILI
			IliDataStruct & , // data
			bool reintegrate=true, // re-integrate
			bool record_lls=false); // store the loglikelihoods of single observations in the data object
	double posterior_loglike(double , // temperature for path sampling
			const double[] , const double[] , double , // background, reportingrate, dispersion
			double , // beta
			double , double , // oscilation of bg ILI
			const double[] , const double[] , // meanLogitS0, sdLogitS0
			double , double , // meanT0, sdT0
			double , double , // meanLogGamma, sdLogGamma
			double , double , // meanLogKappa, sdLogKappa
			double , // sdSeasonalBgEffect
			IliDataStruct & , // data
			bool reintegrate=true); // re-integrate
	IliDataStruct simulate(const double[] , const double[] , double , // background, reporting, dispersion
			double , // beta
			double , double , // oscilation of bg ILI
			const IliDataStruct & , // const data
			Rng & ); // random number generator
	void update(double , // temperature for path sampling
			const double[] , const double[] , double , // background, reportingrate, dispersion
			double , // beta
			double , double , // oscilation of bg ILI
			const double[] , const double[] , // meanLogitS0, sdLogitS0
			double , double , // meanT0, sdT0
			double , double , // meanLogGamma, sdLogGamma
			double , double , // meanLogKappa, sdLogKappa
			double , // sdSeasonalBgEffect
			IliDataStruct & , // data
			Rng & ); // random number generator
	void print(std::ostream & ) const; // for pretty-printing
	std::string xmlString() const;
	int countParameters() const;
protected:
	int seasonNr;
	void updateSinglePar(double , // temperature for path sampling
			const double[] , const double[] , double , // background, reportingrate, dispersion
			double , // beta,
			double , double , // oscilation of bg ILI
			const double[] , const double[] , // meanLogitS0, sdLogitS0
			double , double , // meanT0, sdT0
			double , double , // meanLogGamma, sdLogGamma
			double , double , // meanLogKappa, sdLogKappa
			double , // sdSeasonalBgEffect
			IliDataStruct & , // data
			ParPrior & , Rng & ); // the parameter and a RNG
	void setDefaultValues(); // default values for the parameters / priors / bounds
	// parameters...
	ParPrior S0[NUMBER_OF_AGE_CLASSES];
	ParPrior t0;
	ParPrior gamma; // recovery rate (season specific)
	ParPrior kappa; // incubation rate (season specific)
	ParPrior seasonalBgEffect;
	// integration data..
	double **Ss, **SsA, **SsB; // Ss equals either to SsA or SsB
	int numberUniqueTimePoints; // length of Ss
	int numberOfWeeks; // length of isEpidemicWeek
	std::vector<bool> isEpidemicWeek; // used to determine tStart and tEnd
	// compound parameters
	CompoundPar R0;
	CompoundPar Sinf[NUMBER_OF_AGE_CLASSES]; // computed by integration....
	CompoundPar tStart, tEnd; // used to define an epidemic window
	bool tStartDefined, tEndDefined; // defining tStart and tEnd may be problematic
private:
	void copy(const McmcSeasonState & );
	void clear();
	virtual void mapParameters();
};

std::ostream & operator<<(std::ostream & , const McmcSeasonState & );

// TODO: handle (shared) data better: this is a mess!!!
class SsUpdateJob : public Job {
public:
	SsUpdateJob(McmcSeasonState* , int , // season
		double , // temperature T
		double* , double* , // mean background, slope background
		double* , double* , // slope, offset (reporting)
		double , // dispersion
		double , // beta
		double , double , // ampBg, phaseBg
		double* , double* , double* , // meanLogitS0, sdLogitSo, slopeLogitS0
		double , double , // meanT0, sdT0
		double , double , // meanLogGamma, sdLogGamma
		double , double , // meanLogKappa, sdLogKappa
		double , // sdSeasonalBgEffect
		IliDataStruct* ); // data
protected:
	bool execute(Rng & );
	McmcSeasonState* seasonState;
	int s;
	double T; // temperature for path sampling
	double* meanBg;
	double* slopeBg;
	double* offset;
	double* slope;
	double dispersion;
	double beta;
	double ampBg;
	double phaseBg;
	double* meanLogitS0;
	double* sdLogitS0;
	double* slopeLogitS0;
	double meanT0;
	double sdT0;
	double meanLogGamma;
	double sdLogGamma;
	double meanLogKappa;
	double sdLogKappa;
	double sdSeasonalBgEffect;
	IliDataStruct* data;
};


class SsLoglikeJob : public Job {
public:
	SsLoglikeJob(McmcSeasonState* , int	, // season
		double* , double* , // mean background, slope background
		double* , double* , // offset, slope (reporting)
		double , // dispersion
		double , // beta
		double , double , // ampBg, phaseBg
		bool , bool , // reintegrate, record_lls
		IliDataStruct* ); // data
protected:
	bool execute(Rng & ); // although no RNG is needed here...
	McmcSeasonState* seasonState;
	int s;
	double* meanBg; // mean background,
	double* slopeBg; // slope background
	double* offset; // reporting
	double* slope; // reporting
	double dispersion;
	double beta; // beta
	double ampBg; // amplitude background ILI
	double phaseBg; // pase background ILI
	bool reintegrate;
	bool record_lls;
	IliDataStruct* data;
};

#endif
