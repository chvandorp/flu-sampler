#include "mcmc.hpp"

// methods for Mcmc 

Mcmc::Mcmc(int length, int burnin, int thinning, bool wbic_mode, 
           unsigned long seed, std::string vid) : 
		length(length), burnin(burnin), thinning(thinning), 
		print_thinning(thinning), wbic_mode(wbic_mode), vid(vid) {
	rng.seed(seed);
	workerPool.initWorkerPool(NUMBER_OF_THREADS, rng);
	bic_cache = 0.0; bic_cache_good = false;
	waic_cache = 0.0; waic_cache_good = false;
	wbic_cache = 0.0; wbic_cache_good = false;
	
	beta = 0.0; // todo...
}
Mcmc::~Mcmc() {
	workerPool.waitForWorkerPoolToExit(); 
}

void Mcmc::importData() {
	bool cmat_ok(false), ili_ok(false), timing_ok(false), consult_ok(false), cmats_ok(false);
	beta = 1.0; // TODO
	// todo: handle C as a vector; give gsl access to vector::data
	double* C = new double[NUMBER_OF_AGE_CLASSES * NUMBER_OF_AGE_CLASSES];
	std::string cmatFileName;
	std::string iliFileName;
	std::string consultFileName;
	if ( NUMBER_OF_AGE_CLASSES == 1 ) {
		cmatFileName = "m-matrix-small2-aggr.txt"; // TODO: change to tsv
		iliFileName = "sortedIliAndEstXiaData-aggr.txt";
		// TODO: aggr consultation data
	}
	else {
		cmatFileName = "m-matrix-small2.txt";
		iliFileName = "sortedIliAndEstXiaData.txt";
		consultFileName = "sortedGisData.tsv";
        // cmatFileName = "m-matrix-small-10ageclasses.txt";
		// iliFileName = "sortedIliAndEstXiaData-10ageclasses.txt";
		// consultFileName = "sortedGisData-10ageclasses.tsv";
	}
	std::string timingFileName = "seasons-daynumbers.txt";
	
	cmat_ok = importContactMatrix(DATA_FOLDER + cmatFileName, C);
	if ( !cmat_ok ) throw MsgException("Error in file " + cmatFileName + RIGHT_HERE);
	ili_ok = importIliData(DATA_FOLDER + iliFileName, data);		
	if ( !ili_ok ) throw MsgException("Error in file " + iliFileName + RIGHT_HERE);
	timing_ok = importTimingSeasons(DATA_FOLDER + timingFileName, data);
	if ( !timing_ok ) throw MsgException("Error in file " + timingFileName + RIGHT_HERE);
	// TODO: for the current version, consultation data is not used
	// consult_ok = importGisData(DATA_FOLDER + consultFileName, data);
	// if ( !consult_ok ) throw MsgException("Error in file " + consultFileName + RIGHT_HERE);
	cmats_ok = addContactMatrixToIliData(C, data);
	if ( !cmats_ok ) throw MsgException("Error while creating contact matices" + RIGHT_HERE);
	
	delete[] C;
}

void Mcmc::run() {
	importData();
	
	std::list<Lock> locks;
	bool found_locks = importLocks("parameter-locks.tsv", locks);
	if ( !found_locks ) {
		std::cerr << "WARNING: can't open the lock file" << std::endl;
	}
	
	// create the initial state
	McmcFullState state(data);
	state.setMode(WARMUP_MODE); // redundant: warmup is the default.
	state.setLocks(locks); // select a sub-model...
	state.loglike(beta, data, workerPool);
	
	double T = 1.0;
	if ( wbic_mode ) {
		/** if wbic_mode is true, then we need to sample at a "higher
		* temperature" (of lower T in our parameterization).
		* For computing the WBIC, the parameter T must be set to 1/log(N)
		* where N is the number of observations. See reference
		*    Sumio Watanabe,
		*    A Widely Applicable Bayesian Information Criterion,
		*    Journal of Machine Learning Research 14 (2013) 867-897
		* for more details
		*/
		int N = countIliObservations(data);
		T = 1.0/log(N);
		std::cout << "# WBIC mode: sampling at T = " <<  T << std::endl;
	}
	
	// print the initial state
	std::cout << state << std::endl;
	
	// start the main for loop
	EtaEstimator eta(length + burnin);
	for ( int n = -burnin; n <= length; ++n ) {
		/** TESTING: let the temperature decrease (i.e. T increase from 0 to 1)
		 * to avoid "stuck" parameters.
         * Don't do this in wbic mode, T is already lower.
		 */
		if ( !wbic_mode ) {
			if ( n < 0 ) {
				double x = double(n)/burnin;
				T = 1.0 - x*x;
			}
			else T = 1.0;
		}
		
		if ( n == 0 ) state.setMode(SAMPLING_MODE); // stop changing pvar
		// state.update(beta, gamma, data, rng); // TODO: serial option
		state.update(T, beta, data, rng, workerPool);
		if ( n % thinning == 0 && n >= 0 ) { // save to chain
			state.loglike(beta, data, workerPool, false, true); // record lls of single observations in data objects
			chain.push_back(state); // save the thinned chain after burnin
		}
		if ( n >= 0 ) { // keep track of maximum likelihood parameters
			saveMlState(state);
			/* saves a copy of state in mlState if state.ll_cache > mlState.ml_cache
			 * idem for mplState (for MAP estimate)
			 */
		}
		if ( n % print_thinning == 0 && n != -burnin ) { // for printing to screen
			std::cout << n 
			          << " ETA: " << eta << " " 
			          << "T: " << T << " " 
			          << state << std::endl;
		}
		eta.update(); // update eta each step...
	}
	
	// calculate statistics
	calcBic();
	calcWaic();
	calcWbic();
}

void Mcmc::runSingleSeason(int seasonNr) {
	/** This method takes parameters from the last computed chain,
	 * and then tries to estimate the parameters for an individual 
	 * season.
	 */
	if ( chain.empty() ) {
		throw MsgException("No estimate known for global parameters" + RIGHT_HERE);
	}
    
    std::list<Lock> locks;
	bool found_locks = importLocks("parameter-locks.tsv", locks);
	if ( !found_locks ) {
		std::cerr << "WARNING: can't open the lock file" << std::endl;
	}
	
	// create the initial state
	int numberUniqueTimePoints = data[seasonNr].numberUniqueTimePoints;
	int numberOfWeeks = data[seasonNr].numberOfWeeks;
	McmcSeasonState state(numberOfWeeks, numberUniqueTimePoints, seasonNr);
	state.setMode(WARMUP_MODE); // redundant: warmup is the default.
	state.setLocks(locks); // select a submodel...

	// get a set of global variables
	
	bool ok = true;
	double T = 1;
	double bg[NUMBER_OF_AGE_CLASSES];
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		double slp;
		std::string slpname = "slope_bg_" + std::to_string(a);
		ok = ok && mplState.findValue(slpname, slp);
		double off;
		std::string offname = "bg_" + std::to_string(a);
		ok = ok && mplState.findValue(offname, off);
		bg[a] = calcLinearTrend(off, slp, seasonNr, CENTRAL_SEASON);
	}
	double q[NUMBER_OF_AGE_CLASSES];
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		double slp;
		std::string slpname = "slope_q_" + std::to_string(a);
		ok = ok && mplState.findValue(slpname, slp);
		double off;
		std::string offname = "offset_q_" + std::to_string(a);
		ok = ok && mplState.findValue(offname, off);
		q[a] = calcTrendendProb(off, slp, seasonNr, CENTRAL_SEASON);
	}
	double dispersion;
	ok = ok && mplState.findValue("dispersion", dispersion);
	double ampBg;
	ok = ok && mplState.findValue("amplitude_bg", ampBg); 
	double phaseBg;
	ok = ok && mplState.findValue("phase_bg", phaseBg); 
	double meanLogitS0[NUMBER_OF_AGE_CLASSES];
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		/* 
		double slp;
		std::string slpname = "slope_S0_" + std::to_string(a);
		ok = ok && mplState.findValue(slpname, slp);
		double off;
		std::string offname = "mean_S0_" + std::to_string(a);
		ok = ok && mplState.findValue(offname, off);
		meanLogitS0[a] = calcLinearTrend(off, slp, seasonNr, CENTRAL_SEASON);
		*/
		meanLogitS0[a] = 0.0; // TESTING!!!
	}
	double sdLogitS0[NUMBER_OF_AGE_CLASSES];
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		/*
		std::string parname = "sd_S0_" + std::to_string(a);
		ok = ok && mplState.findValue(parname, sdLogitS0[a]);
		*/
		sdLogitS0[a] = 10.0;
	}
	double meanT0;
	ok = ok && mplState.findValue("mean_t0", meanT0); 
	double sdT0;
	ok = ok && mplState.findValue("sd_t0", sdT0); 
	double meanLogGamma;
	ok = ok && mplState.findValue("mean_log_gamma", meanLogGamma); 
	double sdLogGamma;
	ok = ok && mplState.findValue("sd_log_gamma", sdLogGamma); 
	double meanLogKappa;
	ok = ok && mplState.findValue("mean_log_kappa", meanLogKappa); 
	double sdLogKappa;
	ok = ok && mplState.findValue("sd_log_kappa", sdLogKappa); 
	double sdSeasonalBgEffect;
	ok = ok && mplState.findValue("sd_seasonal_bg", sdSeasonalBgEffect);
	
	if ( !ok ) throw MsgException("could not find a parameter in McmcState" + RIGHT_HERE);
	
	std::cout << "ll = " << state.loglike(bg, q, dispersion, beta, ampBg, phaseBg, data[seasonNr], true) << std::endl;
	std::cout << "lpi = " << state.prior_loglike(meanLogitS0, sdLogitS0, 
												 meanT0, sdT0, 
												 meanLogGamma, sdLogGamma, 
												 meanLogKappa, sdLogKappa, 
												 sdSeasonalBgEffect) << std::endl;
	std::cout << "lpl = " << state.posterior_loglike(T, bg, q, dispersion, beta, ampBg, phaseBg, 
													meanLogitS0, sdLogitS0, 
													meanT0, sdT0, 
													meanLogGamma, sdLogGamma,
													meanLogKappa, sdLogKappa,
													sdSeasonalBgEffect, 
													data[seasonNr], true) << std::endl;	
	// print the initial state
	std::cout << state << std::endl;

	// make a chain
	std::list<McmcSeasonState> ss_chain;
	int ss_length = 2000; // TODO
	int ss_burnin = ss_length; // TODO

	// start the main for loop
	EtaEstimator eta(ss_length + ss_burnin);
	for ( int n = -ss_burnin; n <= ss_length; ++n ) {
		if ( n == 0 ) state.setMode(SAMPLING_MODE); // stop changing pvar
		state.update(T, bg, q, dispersion, beta, ampBg, phaseBg, 
		             meanLogitS0, sdLogitS0, 
		             meanT0, sdT0, 
		             meanLogGamma, sdLogGamma,
		             meanLogKappa, sdLogKappa,
		             sdSeasonalBgEffect, 
		             data[seasonNr], rng);
		if ( n % thinning == 0 && n >= 0 ) { // save to chain
			ss_chain.push_back(state); // save the thinned chain after burnin
		}
		if ( n % print_thinning == 0 && n != -ss_burnin ) { // for printing to screen
			std::cout << n << " ETA: " << eta << " " << state << std::endl;
		}
		eta.update(); // update eta each step...
	}
	
	// write results to files
	
	std::fstream fs; // filestream for ss-chain and ss-sili
	std::stringstream ss; // and a stringstream for the filename
	
	// write the ss_chain to file

	ss << "data/ss-chain-" << vid << "-" << seasonNr << ".xml";
	fs.open(ss.str().c_str(), std::fstream::out);
	fs	<< "<mcmc "
		<< "length='" << ss_length << "' "
		<< "thinning='" << thinning << "' "
		<< "burnin='" << ss_burnin << "' "
		<< "vid='" << vid << "' ";
	fs  << ">" << std::endl;
	
	// chain
	
	fs << "<chain season='" << seasonNr << "' >" << std::endl;
	for ( auto ssit = ss_chain.begin(); ssit != ss_chain.end(); ++ssit ) {
		fs << ssit->xmlString() << std::endl;
	}
	fs << "</chain>" << std::endl;
	fs << "</mcmc>";
	fs.close();
	
	// simulate and print...
	
	ss.str(""); // empties the stringstream
	ss << "data/ss-sili-" << vid << "-" << seasonNr << ".xml";
	fs.open(ss.str().c_str(), std::fstream::out);
	int n = 100; // TODO
	
	fs << "<silis "
	   << "vid='" << vid << "' "
	   << "n='" << n << "' "
	   << ">" << std::endl;
	fs << "<ili >" << std::endl;
	fs << data[seasonNr] << std::endl; // the actual ILI data
	fs << "</ili>" << std::endl;
		
	fs << "<sili >" << std::endl;		
	for ( int i = 0; i < n; ++i ) { // take n samples from the posterior
		// sample a state from the chain
		int idx = rng.integer(ss_chain.size());
		auto cit = ss_chain.begin();
		std::advance(cit, idx);
		// simulate data with the parameters from the sampled state
		auto sdata = cit->simulate(bg, q, dispersion, beta, ampBg, phaseBg, 
								   data[seasonNr], rng);
		fs << sdata << std::endl;
	}
	// closing xml tags
	fs << "</sili>" << std::endl
	   << "</silis>";
	   
	fs.close();
}

void Mcmc::print(std::ostream & os) const {
	os << "<mcmc "
	   << "length='" << length << "' "
	   << "thinning='" << thinning << "' "
	   << "wbic_mode='" << (wbic_mode ? 1 : 0) << "' "
	   << "burnin='" << burnin << "' "
	   << "epsilon='" << INOCULUM_SIZE << "' " // epidemic seeding constant
	   << "num_ageclass='" << NUMBER_OF_AGE_CLASSES << "' "
	   << "bic='" << bic_cache << "' "
	   << "bic_valid='" << (bic_cache_good ? 1 : 0) << "' "
   	   << "waic='" << waic_cache << "' "
	   << "waic_valid='" << (waic_cache_good ? 1 : 0) << "' "
   	   << "wbic='" << wbic_cache << "' "
	   << "wbic_valid='" << (wbic_cache_good && wbic_mode ? 1 : 0) << "' "
	   << "vid='" << vid << "' ";
	os << ">" << std::endl;
	
	// ml state
	os << "<ml_state >" << std::endl
	   << mlState.xmlString() << std::endl
	   << "</ml_state >" << std::endl;

	// mpl state
	os << "<mpl_state >" << std::endl
	   << mplState.xmlString() << std::endl
	   << "</mpl_state >" << std::endl;
	
	os << "<chain >" << std::endl;
	
	for ( auto msit = chain.begin(); msit != chain.end(); ++msit ) {
		os << msit->xmlString() << std::endl;
	}
	os << "</chain>" << std::endl;
	
	os << "</mcmc>";
}

void Mcmc::simulateAndPrint(std::ostream & os, int n) {
	if ( !chain.empty() ) {
		int t = 0; int tmax = 2*n; // for printing progress to screen
		// opening xml tags
		os << "<silis "
		   << "vid='" << vid << "' "
		   << "n='" << n << "' "
		   << ">" << std::endl;
		os << "<ili >" << std::endl;
		os << data << std::endl; // the actual ILI data
		os << "</ili>" << std::endl;
		
		// ml state
		os << "<mlsili >" << std::endl;
		for ( int i = 0; i < n; ++i ) { // simulate with max ll state n times
			std::vector<IliDataStruct> sdata = mlState.simulate(beta, data, rng);
			os << sdata << std::endl;
		}
		os << "</mlsili>" << std::endl;
		
		// mpl state (for MAP estimates)
		os << "<mplsili >" << std::endl;
		for ( int i = 0; i < n; ++i ) { // simulate with max pll state n times
			std::vector<IliDataStruct> sdata = mplState.simulate(beta, data, rng);
			os << sdata << std::endl;
			// keep track of progress
			std::cout << "\rfinished " << t++ << " / " << tmax << "..." << std::flush;
		}
		os << "</mplsili>" << std::endl;

		os << "<sili >" << std::endl;		
		for ( int i = 0; i < n; ++i ) { // take n samples from the posterior
			// sample a state from the chain
			int idx = rng.integer(chain.size());
			std::list<McmcFullState>::iterator cit = chain.begin();
			std::advance(cit, idx);
			// simulate data with the parameters from the sampled state
			std::vector<IliDataStruct> sdata = cit->simulate(beta, data, rng);
			os << sdata << std::endl;
			// keep track of progress
			std::cout << "\r# finished " << t++ << " / " << tmax << "..." << std::flush;
			if ( t == tmax ) std::cout << std::endl;
		}
		// closing xml tags
		os << "</sili>" << std::endl
		   << "</silis>";
	}
	else {
		throw MsgException("can't simulate with an empty chain" + RIGHT_HERE);
	}
	// todo: parallel excecution of simulation
}

double Mcmc::calcBic() {
	std::cout << "# WARNING: in Mcmc::bic the number of parameters is semi-manually set. "
		<< "Please check that countParameters contains all parameters." << std::endl;
	int n = countIliObservations(data);
	std::cout << "# number of observations = " << n << std::endl;
	int k = mlState.countParameters();
	std::cout << "# number of parameters = " << k << std::endl;
	double maxll = 0.0;
	bool cache_ok = mlState.getLoglikeCache(maxll);
	double z = -2*maxll + k * log(n);
	if ( !cache_ok ) {
		std::cerr << "# WARNING: invalid ll_cache during computation of BIC" << std::endl;
	}
	else {
		bic_cache = z;
		bic_cache_good = true;
	}
	return z;
}

double Mcmc::calcWaic() {
	if ( waic_cache_good ) return waic_cache;
	// compute WAIC2 using the lls stored in the data object
	std::list<double> waics;
	for ( auto it = data.begin(); it != data.end(); ++it ) {
		// ILI
		for ( int w = 0; w < it->numberOfWeeks; ++w ) {
			if ( it->weekly_ili[w].age_stratified ) {
				for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
					double elpd = calc_elpd(it->weekly_ili[w].lls_by_age[a]);
					waics.push_back(-2*elpd); // -2 for deviance scale
				}
			}
			else {
				double elpd = calc_elpd(it->weekly_ili[w].lls_total);
				waics.push_back(-2*elpd); // -2 for deviance scale
			}
		}
		// consultation
		if ( it->consultation.observed ) {
			if ( it->consultation.age_stratified ) {
				for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
					double elpd = calc_elpd(it->consultation.lls_by_age[a]);
					waics.push_back(-2*elpd);
				}
			}
			else {
				double elpd = calc_elpd(it->consultation.lls_total);
				waics.push_back(-2*elpd);
			}
		}
	}
	int n = waics.size();
	double z = mean(waics) * n;
	double se = sqrt(n*sample_variance(waics));
	std::cout << "# number of observations: " << n << std::endl;
	std::cout << "# estimated WAIC: " << z << std::endl;
	std::cout << "# standard error WAIC: " << se << std::endl; // todo: save se. part of return?
	// store the computed waic
	waic_cache = z;
	waic_cache_good = true;
	return z; 
}

double Mcmc::calcWbic() {
	if ( wbic_cache_good ) return wbic_cache;
	// compute WBIC using the lls stored in the data object
	std::list<double> wbics;
	for ( auto it = data.begin(); it != data.end(); ++it ) {
		// ILI
		for ( int w = 0; w < it->numberOfWeeks; ++w ) {
			if ( it->weekly_ili[w].age_stratified ) {
				for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
					double ell = mean(it->weekly_ili[w].lls_by_age[a]);
					wbics.push_back(-2*ell); // -2 for deviance scale
				}
			}
			else {
				double ell = mean(it->weekly_ili[w].lls_total);
				wbics.push_back(-2*ell); // -2 for deviance scale
			}
		}
		// consultation
		if ( it->consultation.observed ) {
			if ( it->consultation.age_stratified ) {
				for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
					double ell = mean(it->consultation.lls_by_age[a]);
					wbics.push_back(-2*ell);
				}
			}
			else {
				double ell = mean(it->consultation.lls_total);
				wbics.push_back(-2*ell);
			}
		}
	}
	int n = wbics.size();
	double z = mean(wbics) * n;
	std::cout << "# number of observations: " << n << std::endl;
	std::cout << "# estimated WBIC: " << z << std::endl;
	// store the computed waic
	wbic_cache = z;
	wbic_cache_good = true;
	return z; 
}

void Mcmc::printPointwiseWaic(std::ostream & os) const {
	// TODO: check that this is allowed (exceptions when lls are empty lists??)
	os << "<waic "
	   << "vid='" << vid << "' "
	   << ">" << std::endl;
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		os << data[s].waicXmlString() << std::endl;
	}
	os << "</waic>";
}

void Mcmc::saveMlState(const McmcFullState & state) {
	// maximum likelihood...
	double ll;
	if ( state.getLoglikeCache(ll) ) { // ll passes-by-red
		double maxll; 
		bool maxll_cache_good = mlState.getLoglikeCache(maxll);
		if ( !maxll_cache_good || maxll < ll ) {
			mlState = state; // make a copy of state in mlState
		}
	}
	else {
		std::cerr << "# WARNING: invalid ll cache..." << std::endl;
	}
	// maximum posterior likelihood
	double pll;
	if ( state.getPostLoglikeCache(pll) ) { // pll passes-by-red
		double maxpll;
		bool maxpll_cache_good = mplState.getPostLoglikeCache(maxpll);
		if ( !maxpll_cache_good || maxpll < pll) {
			mplState = state; // make a copy of state in mlState
		}
	}
	else {
		std::cerr << "# WARNING: invalid pll cache..." << std::endl;
	}
}

std::ostream & operator<<(std::ostream & os, const Mcmc & chain) {
	chain.print(os);
	return os;
}

