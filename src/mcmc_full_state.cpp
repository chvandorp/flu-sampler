#include "mcmc_full_state.hpp"

McmcFullState::McmcFullState() {
	setDefaultValues();
	// by default, expect no data.
}

McmcFullState::McmcFullState(const std::vector<IliDataStruct> & data) {
	name = "main_model";
	// get the number of time points from the data
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		int np = data[s].numberUniqueTimePoints;
		int nw = data[s].numberOfWeeks;
		seasonStates[s] = McmcSeasonState(nw, np, s);
	}
	setDefaultValues();
	// the parameter map can now be used (e.g. for finding parameters by name)
	mapParameters();
}

McmcFullState::~McmcFullState() {
	clear();
}
McmcFullState::McmcFullState(const McmcFullState & st) : McmcState(st) {
	copy(st);
}
McmcFullState & McmcFullState::operator=(const McmcFullState & st) {
	if ( this != &st ) {
		clear();
		McmcState::operator=(st); // f*cking C++...
		copy(st);
	}
	return *this;
}
McmcState* McmcFullState::dup() const { // return a pointer to the base class
	return new McmcFullState(*this);
}


void McmcFullState::copy(const McmcFullState & st) {
	/* copy all members manually
	 */
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		seasonStates[s] = st.seasonStates[s];
	}
	// parameters
	meanLogGamma = st.meanLogGamma;
	sdLogGamma = st.sdLogGamma;
	meanLogKappa = st.meanLogKappa;
	sdLogKappa = st.sdLogKappa;

	// age-INdependent parameters
	amplitudeBackground = st.amplitudeBackground;
	phaseBackground = st.phaseBackground;
	dispersion = st.dispersion;
	meanT0 = st.meanT0;
	sdT0 = st.sdT0;
	sdSeasonalBgEffect = st.sdSeasonalBgEffect;
	
	// age-dependent parameters
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		backgrounds[a] = st.backgrounds[a];
		slopeBackgrounds[a] = st.slopeBackgrounds[a];
		slopeReportingRates[a] = st.slopeReportingRates[a];
		offsetReportingRates[a] = st.offsetReportingRates[a];
		meanLogitS0[a] = st.meanLogitS0[a];
		sdLogitS0[a] = st.sdLogitS0[a];
		slopeLogitS0[a] = st.slopeLogitS0[a];
	}
	
	/* it is important that the parameters can be found in 
	 * the parameter map. If mapParameters becomes obsolete,
	 * then perhaps this copy function is obsolete?
	 */
	mapParameters(); // todo... should become deprecated
}

void McmcFullState::clear() {
	/* empty */
}



void McmcFullState::setDefaultValues() {
	/* TODO: 
	 * #-define default values in constants.hpp
	 */
	meanT0 = ParPrior(100.0, 1.0, new NormalPrior(100.0, 1000.0), "mean_t0");
	sdT0 = ParPrior(50.0, 1.0, new HalfNormalPrior(1000.0), "sd_t0");
	sdT0.setLBound(0.0);
	amplitudeBackground = ParPrior(0.0, 0.1, new HalfNormalPrior(100.0), "amplitude_bg");
	amplitudeBackground.setLBound(0.0);
	phaseBackground = ParPrior(0.0, 0.01, new UniformPrior(0.0, 2.0*M_PI), "phase_bg");
	phaseBackground.setBounds(0.0, 2.0*M_PI);
	phaseBackground.setHomotopy(CIRCULAR);
	sdSeasonalBgEffect = ParPrior(1.0, 0.1, new HalfNormalPrior(100.0), "sd_seasonal_bg");
	dispersion = ParPrior(0.01, 0.01, new HalfNormalPrior(100.0), "dispersion");
	dispersion.setLBound(0.0);
	double mlg = -log(DURATION_INFECTIOUS_PERIOD);
	meanLogGamma = ParPrior(mlg, 0.01, new NormalPrior(mlg, 0.1), "mean_log_gamma"); 
	// log(gamma) ~ normal => -log(1/gamma) ~ normal => duration ~ lognormal
	sdLogGamma = ParPrior(0.1, 0.01, new HalfNormalPrior(0.1), "sd_log_gamma");
	sdLogGamma.setLBound(0.0);	
	double mlk = -log(DURATION_EXPOSED_PERIOD);
	meanLogKappa = ParPrior(mlk, 0.01, new NormalPrior(mlk, 0.1), "mean_log_kappa"); // todo: not yet tested
	// log(kappa) ~ normal => -log(1/kappa) ~ normal => duration ~ lognormal
	sdLogKappa = ParPrior(0.1, 0.01, new HalfNormalPrior(0.1), "sd_log_kappa"); 
	sdLogKappa.setLBound(0.0);	
		
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		meanLogitS0[a] = ParPrior(0.0, 1.0, new NormalPrior(0.0, 100.0), "mean_S0_" + std::to_string(a)); 
		slopeLogitS0[a] = ParPrior(0.0, 1.0, new NormalPrior(0.0, 100.0), "slope_S0_" + std::to_string(a));
		sdLogitS0[a] = ParPrior(1.0, 0.1, new HalfNormalPrior(100.0), "sd_S0_" + std::to_string(a));
		sdLogitS0[a].setLBound(0.0);

		backgrounds[a] = ParPrior(-2.0, 0.01, new NormalPrior(0.0, 100.0), "bg_" + std::to_string(a));
		slopeBackgrounds[a] = ParPrior(0.0, 1.0, new NormalPrior(0.0, 100.0), "slope_bg_" + std::to_string(a));
		
		slopeReportingRates[a] = ParPrior(0.0, 1.0, new NormalPrior(0.0, 100.0), "slope_q_" + std::to_string(a));
		offsetReportingRates[a] = ParPrior(0.0, 1.0, new NormalPrior(0.0, 100.0), "offset_q_" + std::to_string(a));
	}
}

void McmcFullState::mapParameters() {
	/* this function should become obsolete in the future,
	 * when parameters are defined/stored in the parameter map by default
	 */
	parameters.clear();
	parameters[meanT0.getName()] = &meanT0;
	parameters[sdT0.getName()] = &sdT0;
	parameters[amplitudeBackground.getName()] = &amplitudeBackground;
	parameters[phaseBackground.getName()] = &phaseBackground;
	parameters[sdSeasonalBgEffect.getName()] = &sdSeasonalBgEffect;
	parameters[dispersion.getName()] = &dispersion;
	parameters[meanLogGamma.getName()] = &meanLogGamma;
	parameters[sdLogGamma.getName()] = &sdLogGamma;
	parameters[meanLogKappa.getName()] = &meanLogKappa;
	parameters[sdLogKappa.getName()] = &sdLogKappa;
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		parameters[meanLogitS0[a].getName()] = &meanLogitS0[a];
		parameters[sdLogitS0[a].getName()] = &sdLogitS0[a];
		parameters[slopeLogitS0[a].getName()] = &slopeLogitS0[a];
		parameters[backgrounds[a].getName()] = &backgrounds[a];
		parameters[slopeBackgrounds[a].getName()] = &slopeBackgrounds[a];
		parameters[slopeReportingRates[a].getName()] = &slopeReportingRates[a];
		parameters[offsetReportingRates[a].getName()] = &offsetReportingRates[a];
	}
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		std::string sname = seasonStates[s].getName();
		submodels[sname] = &seasonStates[s];
	}
}

double McmcFullState::prior_loglike() const {
	double lpi = 0.0;
	lpi += meanT0.loglike();
	lpi += sdT0.loglike();
	lpi += amplitudeBackground.loglike();
	lpi += phaseBackground.loglike();
	lpi += sdSeasonalBgEffect.loglike();
	lpi += dispersion.loglike();
	lpi += meanLogGamma.loglike();
	lpi += sdLogGamma.loglike();
	lpi += meanLogKappa.loglike();
	lpi += sdLogKappa.loglike();

	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		// hyper parameters for S0
		lpi += meanLogitS0[a].loglike();
		lpi += sdLogitS0[a].loglike();
		lpi += slopeLogitS0[a].loglike();
		// background and reporting
		lpi += backgrounds[a].loglike();
		lpi += slopeBackgrounds[a].loglike();
		lpi += slopeReportingRates[a].loglike();
		lpi += offsetReportingRates[a].loglike();
	}
	// get the season priors
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		// compute the season-specific meanLogitS0
		double meanLogitS0vals[NUMBER_OF_AGE_CLASSES];
		double sdLogitS0vals[NUMBER_OF_AGE_CLASSES];
		for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
			meanLogitS0vals[a] = calcLinearTrend(meanLogitS0[a].getValue(), 
					slopeLogitS0[a].getValue(), s, CENTRAL_SEASON);
			sdLogitS0vals[a] = sdLogitS0[a].getValue();
		}
		lpi += seasonStates[s].prior_loglike(meanLogitS0vals, sdLogitS0vals,
				meanT0.getValue(), sdT0.getValue(),
				meanLogGamma.getValue(), sdLogGamma.getValue(),
				meanLogKappa.getValue(), sdLogKappa.getValue(),
				sdSeasonalBgEffect.getValue());
	}
	return lpi;
}

// serial method
double McmcFullState::loglike(double beta, 
		std::vector<IliDataStruct> & data, bool reintegrate, bool record_lls) {
	std::cerr << "# WARNING: McmcFullState::loglike (serial) called..." << std::endl; // TODO
	double ll = 0.0;
	// seasons...
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) { // todo: multiprocessing
		double q[NUMBER_OF_AGE_CLASSES];
		double bg[NUMBER_OF_AGE_CLASSES];
		for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
			q[a] = calcTrendendProb(offsetReportingRates[a].getValue(), 
					slopeReportingRates[a].getValue(), s, CENTRAL_SEASON);
			bg[a] = calcLinearTrend(backgrounds[a].getValue(),
			        slopeBackgrounds[a].getValue(), s, CENTRAL_SEASON);
		}
		double sll = seasonStates[s].loglike(bg, q, dispersion.getValue(), beta,
				amplitudeBackground.getValue(), phaseBackground.getValue(),
				data[s], reintegrate, record_lls);
		ll += sll;
	}
	setLoglikeCache(ll);
	return ll;
}

// parallel method
double McmcFullState::loglike(double beta, 
		std::vector<IliDataStruct> & data, WorkerPool & workerPool,
		bool reintegrate, bool record_lls) {
	// get values of parameters that hold for all seasons
	double meanBg[NUMBER_OF_AGE_CLASSES];
	double slopeBg[NUMBER_OF_AGE_CLASSES];
	double offset[NUMBER_OF_AGE_CLASSES];
	double slope[NUMBER_OF_AGE_CLASSES];
	
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		// background and reporting
		meanBg[a] = backgrounds[a].getValue();
		slopeBg[a] = slopeBackgrounds[a].getValue();
		offset[a] = offsetReportingRates[a].getValue();
		slope[a] = slopeReportingRates[a].getValue();		
	}		
			
	// feed jobs to the workerpool
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		// make a new job
		SsLoglikeJob* job = new SsLoglikeJob(&seasonStates[s], s,
		meanBg,	slopeBg, offset, slope,	dispersion.getValue(), beta,
		amplitudeBackground.getValue(), phaseBackground.getValue(), 
		reintegrate, record_lls, &data[s]);
		// job should be deleted when finished
		job->setDeleteWhenFinished();
		workerPool.addNewJob(job);
	}
	// wait for all Jobs to be finished
	workerPool.syncWorkerThreads();

	double ll = 0.0;
	// get season-lls from the caches
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		double sll = 0.0;
		bool ok = seasonStates[s].getLoglikeCache(sll);
		if ( !ok ) {
			std::cerr << "WARNING: loglike cache not stored correctly." << std::endl;
		}
		else {
			ll += sll;
		}
	}
	setLoglikeCache(ll);
	return ll;
}

// parallel method
double McmcFullState::posterior_loglike(double T,
		double beta, 
		std::vector<IliDataStruct> & data, WorkerPool & workerPool,
		bool reintegrate) {
	double pll = prior_loglike();
	if ( T != 0 ) {
		pll += T * loglike(beta, data, workerPool, reintegrate);
	}
	setPostLoglikeCache(pll);
	return pll;
}

// serial method
double McmcFullState::posterior_loglike(double T,
		double beta, 
		std::vector<IliDataStruct> & data, bool reintegrate) {
	double pll = prior_loglike();
	if ( T != 0 ) {
		pll += T * loglike(beta, data, reintegrate);
	}
	setPostLoglikeCache(pll);
	return pll;
}


std::vector<IliDataStruct> McmcFullState::simulate(double beta,
		const std::vector<IliDataStruct> & data,
		Rng & rng) {
	std::vector<IliDataStruct> sdata;
	sdata.reserve(NUMBER_OF_SEASONS);
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		// get values of parameters that are season specific
		double q[NUMBER_OF_AGE_CLASSES];
		double bg[NUMBER_OF_AGE_CLASSES];
		for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
			q[a] = calcTrendendProb(offsetReportingRates[a].getValue(), 
					slopeReportingRates[a].getValue(), s, CENTRAL_SEASON);
			bg[a] = calcLinearTrend(backgrounds[a].getValue(),
			        slopeBackgrounds[a].getValue(), s, CENTRAL_SEASON);
		}
		IliDataStruct ssdata = seasonStates[s].simulate(bg, q, 
				dispersion.getValue(), beta,
				amplitudeBackground.getValue(), phaseBackground.getValue(),
				data[s], rng);
		sdata.push_back(ssdata);
	}
	return sdata;
}


// serial update method
void McmcFullState::update(double T, double beta, 
		std::vector<IliDataStruct> & data, Rng & rng) {
	std::cerr << "# WARNING: McmcFullState::update (serial) has not been tested properly..." << std::endl; // TODO
	// update global parameters...
	updateSinglePar(T, beta, data, meanT0, rng, false);
	updateSinglePar(T, beta, data, sdT0, rng, false);
	updateSinglePar(T, beta, data, amplitudeBackground, rng, false);
	updateSinglePar(T, beta, data, phaseBackground, rng, false);
	updateSinglePar(T, beta, data, sdSeasonalBgEffect, rng, false);
	updateSinglePar(T, beta, data, dispersion, rng, false);
	updateSinglePar(T, beta, data, meanLogGamma, rng, false);
	updateSinglePar(T, beta, data, sdLogGamma, rng, false);
	updateSinglePar(T, beta, data, meanLogKappa, rng, false);
	updateSinglePar(T, beta, data, sdLogKappa, rng, false);
	
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		// hyper parameters for S0
		updateSinglePar(T, beta, data, meanLogitS0[a], rng, false);
		updateSinglePar(T, beta, data, sdLogitS0[a], rng, false);
		updateSinglePar(T, beta, data, slopeLogitS0[a], rng, false);
		// background and reporting
		updateSinglePar(T, beta, data, backgrounds[a], rng, false);
		updateSinglePar(T, beta, data, slopeBackgrounds[a], rng, false);
		updateSinglePar(T, beta, data, slopeReportingRates[a], rng, false);
		updateSinglePar(T, beta, data, offsetReportingRates[a], rng, false);
	}
	
	// update season-specific parameters...
	
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		// get values of parameters that are season specific
		double q[NUMBER_OF_AGE_CLASSES];
		double bg[NUMBER_OF_AGE_CLASSES];
		double meanLogitS0vals[NUMBER_OF_AGE_CLASSES];
		double sdLogitS0vals[NUMBER_OF_AGE_CLASSES];		
		for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
			// reporting
			q[a] = calcTrendendProb(offsetReportingRates[a].getValue(), 
					slopeReportingRates[a].getValue(), s, CENTRAL_SEASON);
			// background ILI
			bg[a] = calcLinearTrend(backgrounds[a].getValue(),
			        slopeBackgrounds[a].getValue(), s, CENTRAL_SEASON);
			// susceptiblity (mixed effects model)
			meanLogitS0vals[a] = calcLinearTrend(meanLogitS0[a].getValue(),
					slopeLogitS0[a].getValue(), s, CENTRAL_SEASON);
			sdLogitS0vals[a] = sdLogitS0[a].getValue();
		}
		seasonStates[s].update(T, bg, q, dispersion.getValue(), beta,
				amplitudeBackground.getValue(), phaseBackground.getValue(),
				meanLogitS0vals, sdLogitS0vals,
				meanT0.getValue(), sdT0.getValue(),
				meanLogGamma.getValue(), sdLogGamma.getValue(), 
				meanLogKappa.getValue(), sdLogKappa.getValue(), 
				sdSeasonalBgEffect.getValue(), 
				data[s], rng);
	}
}

// parallel method
void McmcFullState::update(double T, 
		double beta,
		std::vector<IliDataStruct> & data, 
		Rng & rng, WorkerPool & workerPool) {
	// update global parameters...
	updateSinglePar(T, beta, data, meanT0, rng, workerPool, false);
	updateSinglePar(T, beta, data, sdT0, rng, workerPool, false);
	updateSinglePar(T, beta, data, amplitudeBackground, rng, workerPool, false);
	updateSinglePar(T, beta, data, phaseBackground, rng, workerPool, false);
	updateSinglePar(T, beta, data, sdSeasonalBgEffect, rng, workerPool, false);
	updateSinglePar(T, beta, data, dispersion, rng, workerPool, false);
	updateSinglePar(T, beta, data, meanLogGamma, rng, workerPool, false);
	updateSinglePar(T, beta, data, sdLogGamma, rng, workerPool, false);
	updateSinglePar(T, beta, data, meanLogKappa, rng, workerPool, false);
	updateSinglePar(T, beta, data, sdLogKappa, rng, workerPool, false);


	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		// S0 hyper parameters
		updateSinglePar(T, beta, data, meanLogitS0[a], rng, workerPool, false);
		updateSinglePar(T, beta, data, sdLogitS0[a], rng, workerPool, false);
		updateSinglePar(T, beta, data, slopeLogitS0[a], rng, workerPool, false);
		// background and reporting
		updateSinglePar(T, beta, data, backgrounds[a], rng, workerPool, false);
		updateSinglePar(T, beta, data, slopeBackgrounds[a], rng, workerPool, false);
		updateSinglePar(T, beta, data, slopeReportingRates[a], rng, workerPool, false);
		updateSinglePar(T, beta, data, offsetReportingRates[a], rng, workerPool, false);
	}
	
	// update season-specific parameters...
	
	// get values of parameters that hold for all seasons
	double meanBg[NUMBER_OF_AGE_CLASSES];
	double slopeBg[NUMBER_OF_AGE_CLASSES];
	double offset[NUMBER_OF_AGE_CLASSES];
	double slope[NUMBER_OF_AGE_CLASSES];
	
	double meanLogitS0vals[NUMBER_OF_AGE_CLASSES];
	double sdLogitS0vals[NUMBER_OF_AGE_CLASSES];
	double slopeLogitS0vals[NUMBER_OF_AGE_CLASSES];
	
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		// the season-specific S0 priors have to be computed by the Job...
		meanLogitS0vals[a] = meanLogitS0[a].getValue();
		sdLogitS0vals[a] = sdLogitS0[a].getValue();
		slopeLogitS0vals[a] = slopeLogitS0[a].getValue();
		// background and reporting
		meanBg[a] = backgrounds[a].getValue();
		slopeBg[a] = slopeBackgrounds[a].getValue();
		offset[a] = offsetReportingRates[a].getValue();
		slope[a] = slopeReportingRates[a].getValue();		
	}
	// create a Job for each season
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		SsUpdateJob* job = new SsUpdateJob(&seasonStates[s], s, // NEW
			T, meanBg, slopeBg, offset, slope, dispersion.getValue(), beta,
			amplitudeBackground.getValue(), phaseBackground.getValue(),
			meanLogitS0vals, sdLogitS0vals, slopeLogitS0vals, 
			meanT0.getValue(), sdT0.getValue(),
			meanLogGamma.getValue(), sdLogGamma.getValue(),
			meanLogKappa.getValue(), sdLogKappa.getValue(),
			sdSeasonalBgEffect.getValue(), 
			&data[s]);
		job->setDeleteWhenFinished(); // sets deleteWhenFinished to true
		workerPool.addNewJob(job);
	}
	workerPool.syncWorkerThreads();
	/* WorkerPool deletes jobs when they are finished 
	 * (if deleteWhenFinished flag is set to true)
	 */
}

// serial method
void McmcFullState::updateSinglePar(double T, double beta,
		std::vector<IliDataStruct> & data, 
		ParPrior & par, Rng & rng, bool reintegrate) {
	double pll_old = posterior_loglike(T, beta, data, false); // makes use of the old Ss
	par.mutate(rng);
	double pll_new = posterior_loglike(T, beta, data, reintegrate);
	double logu = log(rng.uniform());
	if ( pll_new >= pll_old + logu ) {
		par.ac1up(); // accept!
	}
	else {
		par.go_back(); // reject!
		resetLoglikeCache();
		resetPostLoglikeCache();
	}
	if ( mode == WARMUP_MODE ) {
		par.updatePvar();
	}
}

// parallel method
void McmcFullState::updateSinglePar(double T, double beta,
		std::vector<IliDataStruct> & data, 
		ParPrior & par, Rng & rng, WorkerPool & workerPool, 
		bool reintegrate) {
	double pll_old = posterior_loglike(T, beta, data, workerPool, false); // makes use of the old Ss
	par.mutate(rng);
	double pll_new = posterior_loglike(T, beta, data, workerPool, reintegrate);
	double logu = log(rng.uniform());
	if ( pll_new >= pll_old + logu ) {
		par.ac1up(); // accept!
	}
	else {
		par.go_back(); // reject!
		resetLoglikeCache();
		resetPostLoglikeCache();
	}
	if ( mode == WARMUP_MODE ) {
		par.updatePvar();
	}
}

void McmcFullState::print(std::ostream & os) const {
	std::string red = BRED_TERM_FONT;
	std::string def = DEF_TERM_FONT;
	
	// store current precision, and restore before returning
	std::streamsize stream_prec = os.precision(); 
	os.precision(3);

	os << ( !ll_cache_good ? red : "" ) << "ll: " << ll_cache << " " << def << std::endl;
	os << "S0 hyper params (mean):";
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		os << " " << meanLogitS0[a]; 
	}
	os << std::endl;
	os << "S0 hyper params (sd):";
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		os << " " << sdLogitS0[a]; 
	}
	os << std::endl;
	os << "S0 hyper params (slope):";
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		os << " " << slopeLogitS0[a]; 
	}
	os << std::endl;
	os << "t0 hyper params: mean: " << meanT0 << " sd: " << sdT0 << std::endl;
	os << "gamma hyper params: mean: " << meanLogGamma << " sd: " << sdLogGamma << std::endl;
	os << "kappa hyper params: mean: " << meanLogKappa << " sd: " << sdLogKappa << std::endl;
	os << "dispersion: " << dispersion << std::endl;
	os << "bg ili: age effects: ";
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		os << " " << backgrounds[a];
	}
	os << std::endl;
	os << "bg ili: slopes: ";
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		os << " " << slopeBackgrounds[a];
	}
	os << std::endl;
	os << "bg ili: " 
	   << "amplitude: " << amplitudeBackground << " " 
	   << "phase: " << phaseBackground << " " 
	   << "sd: " << sdSeasonalBgEffect << std::endl;
	os << "reportingrate (slope): ";
		for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		os << " " << slopeReportingRates[a];
	}
	os << std::endl;
	os << "reportingrate (offset): "; 
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		os << " " << offsetReportingRates[a];
	}
	os << std::endl;
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		os << s << " " << seasonStates[s] << std::endl;
	}
	
	os.precision(stream_prec); // restore...
}

std::string McmcFullState::xmlString() const {
	std::stringstream ss;
	ss << "<mcmc_state "
	   << "ll='" << ll_cache << "' "
	   << "ll_valid='" << (ll_cache_good ? 1 : 0) << "' "
   	   << "pll='" << pll_cache << "' "
	   << "pll_valid='" << (pll_cache_good ? 1 : 0) << "' "
	   << "k='" << countParameters() << "' "
	   << ">" << std::endl; // end of opening tag 'mcmc_state'
	// print parameters' xml strings
	ss << meanT0.xmlString() << std::endl;
	ss << sdT0.xmlString() << std::endl;
	ss << amplitudeBackground.xmlString() << std::endl;
	ss << phaseBackground.xmlString() << std::endl;
	ss << sdSeasonalBgEffect.xmlString() << std::endl;
	ss << dispersion.xmlString() << std::endl;
	ss << meanLogGamma.xmlString() << std::endl;
	ss << sdLogGamma.xmlString() << std::endl;
	ss << meanLogKappa.xmlString() << std::endl;
	ss << sdLogKappa.xmlString() << std::endl;
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		ss << meanLogitS0[a].xmlString() << std::endl;
		ss << sdLogitS0[a].xmlString() << std::endl;
		ss << slopeLogitS0[a].xmlString() << std::endl;
		ss << backgrounds[a].xmlString() << std::endl;
		ss << slopeBackgrounds[a].xmlString() << std::endl;
		ss << slopeReportingRates[a].xmlString() << std::endl;
		ss << offsetReportingRates[a].xmlString() << std::endl;
	}
	// season-specific parameters...
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		ss << seasonStates[s].xmlString() << std::endl;
	}
	ss << "</mcmc_state>";
	return ss.str();
}

std::ostream & operator<<(std::ostream & os, const McmcFullState & state) {
	state.print(os);
	return os;
}

int McmcFullState::countParameters() const {
	// TODO/FIXME: use the parameter map to count parameters...
	int k = 0;
	
	// single globals	
	if ( !amplitudeBackground.isLocked() ) ++k;
	if ( !phaseBackground.isLocked() ) ++k;
	if ( !sdSeasonalBgEffect.isLocked() ) ++k;
	if ( !dispersion.isLocked() ) ++k;
	if ( !meanT0.isLocked() ) ++k;
	if ( !sdT0.isLocked() ) ++k;
	if ( !meanLogGamma.isLocked() ) ++k;
	if ( !sdLogGamma.isLocked() ) ++k;
	if ( !meanLogKappa.isLocked() ) ++k;
	if ( !sdLogKappa.isLocked() ) ++k;

	// age-dependent globals
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		if ( !meanLogitS0[a].isLocked() ) ++k;
		if ( !sdLogitS0[a].isLocked() ) ++k;
		if ( !slopeLogitS0[a].isLocked() ) ++k;

		if ( !backgrounds[a].isLocked() ) ++k;
		if ( !slopeBackgrounds[a].isLocked() ) ++k;
		if ( !offsetReportingRates[a].isLocked() ) ++k;
		if ( !slopeReportingRates[a].isLocked() ) ++k;
	}
	
	// season-dependent (count parameters of sub models)
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		k += seasonStates[s].countParameters();
	}
	
	return k;
}
