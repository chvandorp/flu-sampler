#include "mcmc_season_state.hpp"

McmcSeasonState::McmcSeasonState() {
	// set default values for parameters / priors / bounds
	setDefaultValues();
	Ss = NULL; // expect no data by default
	SsA = NULL; SsB = NULL;
	numberUniqueTimePoints = 0;
	numberOfWeeks = 0;
	seasonNr = 0;
	tStartDefined = false;
	tEndDefined = false;
}
McmcSeasonState::McmcSeasonState(int numberOfWeeks, int numberUniqueTimePoints, int seasonNr) {
	this->seasonNr = seasonNr;
	name = "sub_model_" + std::to_string(seasonNr);
	// set default values for parameters / priors / bounds
	setDefaultValues();
	// the parameter map can be used.. careful with copying...
	mapParameters();
	// init integration data array
	this->numberUniqueTimePoints = numberUniqueTimePoints;
	this->numberOfWeeks = numberOfWeeks;
	SsA = new double*[numberUniqueTimePoints];
	SsB = new double*[numberUniqueTimePoints];
	for ( int w = 0; w < numberUniqueTimePoints; ++w ) {
		SsA[w] = new double[NUMBER_OF_AGE_CLASSES];
		SsB[w] = new double[NUMBER_OF_AGE_CLASSES];
		for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
			SsA[w][a] = 1.0; // no epidemic by default.
			SsB[w][a] = 1.0;
		}
	}
	Ss = SsA; // let Ss point at one of the caches
	// initiate isEpidemicWeek used to determine tStart and tStop
	isEpidemicWeek.resize(numberOfWeeks, false);
	tStartDefined = false;
	tEndDefined = false;
}
McmcSeasonState::~McmcSeasonState() {
	clear();
}
McmcSeasonState::McmcSeasonState(const McmcSeasonState & sist) : McmcState(sist) {
	copy(sist);
}
McmcSeasonState & McmcSeasonState::operator=(const McmcSeasonState & sist) {
	if ( this != &sist ) {
		clear();
		McmcState::operator=(sist); // f*cking C++...
		copy(sist);
	}
	return *this;
}
McmcState* McmcSeasonState::dup() const { // return a pointer to the base class
	return new McmcSeasonState(*this);
}

void McmcSeasonState::copy(const McmcSeasonState & sist) {
	// copy all members manually
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		S0[a] = sist.S0[a];
	}
	t0 = sist.t0;
	gamma = sist.gamma;
	kappa = sist.kappa;
	seasonalBgEffect = sist.seasonalBgEffect;
	// copy caches
	numberUniqueTimePoints = sist.numberUniqueTimePoints;
	numberOfWeeks = sist.numberOfWeeks;
	seasonNr = sist.seasonNr;
	SsA = new double*[numberUniqueTimePoints];
	SsB = new double*[numberUniqueTimePoints];
	for ( int w = 0; w < numberUniqueTimePoints; ++w ) {
		SsA[w] = new double[NUMBER_OF_AGE_CLASSES];
		SsB[w] = new double[NUMBER_OF_AGE_CLASSES];
		for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
			SsA[w][a] = sist.Ss[w][a];
			SsB[w][a] = 1.0; // just default the other cache to disease free
		}
	}
	Ss = SsA;
	// copy computed values
	R0 = sist.R0;
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		Sinf[a] = sist.Sinf[a];
	}
	tStart = sist.tStart;
	tEnd = sist.tEnd;
	isEpidemicWeek = sist.isEpidemicWeek;
	tStartDefined = sist.tStartDefined;
	tEndDefined = sist.tEndDefined;
	mapParameters(); // todo... should become obsolete
}
void McmcSeasonState::clear() {
	if ( SsA != NULL ) {
		for ( int w = 0; w < numberUniqueTimePoints; ++w ) {
			delete[] SsA[w];
		}
		delete[] SsA;
		SsA = NULL;
	}
	if ( SsB != NULL ) {
		for ( int w = 0; w < numberUniqueTimePoints; ++w ) {
			delete[] SsB[w];
		}
		delete[] SsB;
		SsB = NULL;
	}	
	Ss = NULL;
}

void McmcSeasonState::setDefaultValues() {
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		std::stringstream S0sstr;
		S0sstr << "S0_" << seasonNr << "_" << a;
		S0[a] = ParPrior(0.1, 0.001, new StdNormalPrior, S0sstr.str()); // hypo parameter
		S0[a].setBounds(0.0, 1.0);
	}
	std::string t0str = "t0_" + std::to_string(seasonNr);
	t0 = ParPrior(100.0, 1.0, new StdNormalPrior, t0str); // hypo parameter
	t0.setBounds(0.0, DAYS_IN_YEAR); // fixme: prior has support outside parameter domain
	
	std::string gammastr = "gamma_" + std::to_string(seasonNr);
	gamma = ParPrior(1/DURATION_INFECTIOUS_PERIOD, 0.01, new StdNormalPrior, gammastr); // hypo parameter
	gamma.setLBound(0.0);
	
	std::string kappastr = "kappa_" + std::to_string(seasonNr);
	kappa = ParPrior(1/DURATION_EXPOSED_PERIOD, 0.01, new StdNormalPrior, kappastr); // hypo parameter
	kappa.setLBound(0.0);

	std::string bgstr = "seasonal_bg_" + std::to_string(seasonNr);
	seasonalBgEffect = ParPrior(0, 0.1, new StdNormalPrior, bgstr); // hypo parameter.
	// todo: #-define default values in constants.hpp
	
	// set computed values
	std::string R0str = "R0_" + std::to_string(seasonNr);
	R0 = CompoundPar(0.0, R0str);
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		std::stringstream Sinfstr;
		Sinfstr << "Sinf_" << seasonNr << "_" << a;
		Sinf[a] = CompoundPar(S0[a].getValue(), Sinfstr.str());
	}
	std::string tStartStr = "t_start_" + std::to_string(seasonNr);
	tStart = CompoundPar(0.0, tStartStr);
	std::string tEndStr = "t_end_" + std::to_string(seasonNr);
	tEnd = CompoundPar(0.0, tEndStr);
}

void McmcSeasonState::mapParameters() {
	parameters.clear();
	parameters[t0.getName()] = &t0;
	parameters[gamma.getName()] = &gamma;
	parameters[kappa.getName()] = &kappa;
	parameters[seasonalBgEffect.getName()] = &seasonalBgEffect;
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		parameters[S0[a].getName()] = &S0[a];
	}
}
double McmcSeasonState::prior_loglike(const double meanLogitS0[], const double sdLogitS0[],
		double meanT0, double sdT0, 
		double meanLogGamma, double sdLogGamma, 
		double meanLogKappa, double sdLogKappa, 
		double sdSeasonalBgEffect) const {
	double lpi = 0.0;
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		LogitTransformation fun(meanLogitS0[a], sdLogitS0[a]); // S0 has a 'hypo-prior'
		lpi += S0[a].loglike(&fun);
	}
	LinearTransformation gun(meanT0, sdT0); // t0 is a 'hypo-prior'
	lpi += t0.loglike(&gun);
	
	LinearTransformation hun(0.0, sdSeasonalBgEffect); // soft centering
	lpi += seasonalBgEffect.loglike(&hun); // todo: hyper parameter...
	/* FIXME: soft-centering seasonal background does not work.
	 * too much correlated with age-dependent background
	 */

	LogTransformation logtransform_gamma(meanLogGamma, sdLogGamma); // gamma has a 'hypo-prior'
	lpi += gamma.loglike(&logtransform_gamma);

	LogTransformation logtransform_kappa(meanLogKappa, sdLogKappa); // kappa has a 'hypo-prior'
	lpi += kappa.loglike(&logtransform_kappa);

	return lpi;
}
double McmcSeasonState::loglike(const double bg[], const double q[], // noise and reporting
		double dispersion, // dispersion (poisson or polya ILI cases?)
		double beta, // epidemic parameters
		double ampBg, double phaseBg, // oscilation of bg ILI
		IliDataStruct & data, // data
		bool reintegrate, // if reintegrate is true, do the integration again...
		bool record_lls) { // if record_lls is true, push back ll_obs in data object
	double ll = 0.0;
	
	OdeParamStruct opars;
	opars.C = data.C;
	opars.S0 = new double[NUMBER_OF_AGE_CLASSES];
	opars.beta = beta;
	opars.gamma = gamma.getValue();
	opars.kappa = kappa.getValue();
	opars.t0 = t0.getValue() + data.firstDay;
	opars.ni = NUMBER_OF_INFECTIOUS_STAGES; 
	opars.ne = NUMBER_OF_EXPOSED_STAGES; 
	opars.na = NUMBER_OF_AGE_CLASSES;
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		opars.S0[a] = S0[a].getValue();
	}
	
	double* ts = data.ts;
	int len_ts = data.numberUniqueTimePoints; // fixme: the class has this as a member
	int weeks = data.numberOfWeeks; // fixme: idem
	
	R0 = calcR0(opars); 
	
	if ( reintegrate ) {
		if ( R0 > 1.0 ) {
			// integrate_model2siir(ts, Ss, len_ts, pars);
            // integrate_model2seir(ts, Ss, len_ts, pars);
            integrate_general_model(ts, Ss, len_ts, opars);
		}
		else { 
			// just fill Ss with S0 (warning: we get a minor discontinuity in the likelihood)
			for ( int i = 0; i < len_ts; ++i ) {
				for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
					Ss[i][a] = S0[a].getValue();
				}
			}
		}
		for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
			/* Sinf is used to compute attack rate,
			 * without solving a difficult final size equation
			 * NB: the result may deviate from the "real" final size
			 * FIXME: sometimes we still get negative attack rates...
			 */
			Sinf[a] = Ss[numberUniqueTimePoints-1][a];
		}
	}
		
	int i = 0; // counter to run through Ss and ts
	for ( int w = 0; w < weeks; ++w, ++i ) {
		if ( i >= len_ts ) {
			std::cerr << "WARNING. about to run out of time points..." << std::endl;
			break;
		}
		double totalIncidence = 0.0; // i.e. flu
		int totalPopsize = 0; // sum of g_1, g_2, ...
		double firstDay = data.weekly_ili[w].firstDay;
		double lastDay = data.weekly_ili[w].lastDay;
		bool age_stratified = data.weekly_ili[w].age_stratified;
		int f_tot = 0; // used for non-age-stratified data
		double lambda_tot = 0.0; // used for non-age-stratified data
		if ( ts[i] != firstDay ) ++i; // FIXME: ????
		if ( ts[i] != firstDay || ts[i+1] != lastDay ) {
			std::cerr << "WARNING. time mismatch! " << i << ": " 
					<< ts[i] << " /= " << firstDay << " or "
					<< ts[i+1] << " /= " << lastDay << std::endl;
		}
		for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
			int f = data.weekly_ili[w].gf_by_age[a].second; // (g,f)
			int g = data.weekly_ili[w].gf_by_age[a].first;
			double DeltaS = std::max(0.0, Ss[i][a] - Ss[i+1][a]); // 0.0 if R0 < 1
			double bg_se = seasonalBgEffect.getValue();
			double bgili = calcBackgroundIli(bg[a]+bg_se, ampBg, phaseBg, firstDay);
			double lambda = g * q[a] * (DeltaS + bgili - DeltaS*bgili);
			if ( age_stratified ) {
				double ll_obs = ll_polya(lambda, dispersion, f);
				ll += ll_obs;
				if ( record_lls ) {
					data.weekly_ili[w].lls_by_age[a].push_back(ll_obs); // for WAIC
				}
			}
			else {
				lambda_tot += lambda;
				f_tot += f;
			}
			totalIncidence += DeltaS * g;
			totalPopsize += g;
		}
		if ( !age_stratified ) {
			if ( dispersion > 0.0 ) {
				std::cerr << "WARNING: NegBin distributed ran vars are not additive" 
				          << std::endl;
				// FIXME...
			}
			double ll_obs = ll_polya(lambda_tot, dispersion, f_tot);
			ll += ll_obs;
			if ( record_lls ) {
				data.weekly_ili[w].lls_total.push_back(ll_obs); // for WAIC
			}
		}
		// is this week part of the epidemic?
		if ( totalIncidence / totalPopsize >= EPIDEMIC_THRESHOLD ) {
			isEpidemicWeek[w] = true;
		}
		else { 
			isEpidemicWeek[w] = false;
		}
	}
	
	/* use the boolean vector isEpidemicWeek to determine tStart and tEnd
	 * isEpidemicWeek ideally has the following form:
	 * 00000000111110000
	 * but the following might also happen:
	 * 000001011001111101000
	 * tStart is the first day of the first epidemic week,
	 * tEnd is the first day of the first non-epidemic week
	 */
	tStart = 0.0;
	tEnd = 0.0;
	tStartDefined = false;
	tEndDefined = false;
	for ( int w = 0; w < weeks; ++w ) {
		if ( !tStartDefined && isEpidemicWeek[w] ) {
			tStart = data.weekly_ili[w].firstDay;
			tStartDefined = true;
		}
		if ( !tEndDefined && !isEpidemicWeek[w] && tStartDefined ) {
			tEnd = data.weekly_ili[w].firstDay;
			tEndDefined = true;
		}
	}
	
	/* compute loglikelihood of consulation data.
	 * i.e. the number of GP visits d amoung ILI+ patients f
	 * is distributed d ~ Binomial(f,q)
	 */
	if ( data.consultation.observed && data.consultation.age_stratified ) {
		for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
			int f = data.consultation.fd_by_age[a].first;
			int d = data.consultation.fd_by_age[a].second;
			double ll_obs = ll_binomial(f, q[a], d);
			ll += ll_obs;
			if ( record_lls ) { // for WAIC
				data.consultation.lls_by_age[a].push_back(ll_obs);
			}
		}
	}
	
	// clean up
	delete[] opars.S0;
	// save the ll to ll_cache (saves the previous value...)
	setLoglikeCache(ll);
	return ll;
}

double McmcSeasonState::posterior_loglike(double T, // temperature for path sampling
		const double bg[], const double q[], // noise and reporting
		double dispersion,
		double beta, // epidemic parameters
		double ampBg, double phaseBg, // oscilation of bg ILI
		const double meanLogitS0[], const double sdLogitS0[], // hyper parameters
		double meanT0, double sdT0, // hyper parameters
		double meanLogGamma, double sdLogGamma, // recovery rate hypers
		double meanLogKappa, double sdLogKappa, // incubation rate hypers
		double sdSeasonalBgEffect, // background ILI hypers
		IliDataStruct & data, // data
		bool reintegrate) { // if true, do the integration again...
	double pll = prior_loglike(meanLogitS0, sdLogitS0, meanT0, sdT0, 
						meanLogGamma, sdLogGamma, meanLogKappa, sdLogKappa, 
						sdSeasonalBgEffect);
	if ( T != 0 ) {
		pll += T * loglike(bg, q, dispersion, beta, ampBg, phaseBg, data, reintegrate);
	}
	setPostLoglikeCache(pll);
	return pll;
}


IliDataStruct McmcSeasonState::simulate(const double bg[], const double q[],
		double dispersion,
		double beta,
		double ampBg, double phaseBg, // oscilation of bg ILI
		const IliDataStruct & data,
		Rng & rng) {
	IliDataStruct sdata = data; // TODO: warning! no proper copy constructor;
	sdata.clear_lls(); // TODO/FIXME: this is just a bodge
	sdata.ts = NULL; 
	sdata.C = NULL; // TODO/FIXME: automate with a copy constructor(?)
	// not making ts an C NULL will give a segfault

	double* ts = data.ts;
	int len_ts = data.numberUniqueTimePoints;
	int weeks = data.numberOfWeeks;
	
	// ILI
	
	int i = 0; // counter to run through Ss and ts
	for ( int w = 0; w < weeks; ++w, ++i ) { // lots of duplicate code.. todo (see loglike)
		if ( i >= len_ts ) {
			std::cerr << "about to run out of time points..." << std::endl;
			break;
		}
		double firstDay = data.weekly_ili[w].firstDay;
		double lastDay = data.weekly_ili[w].lastDay;
		if ( ts[i] != firstDay ) ++i;
		if ( ts[i] != firstDay || ts[i+1] != lastDay ) {
			std::cerr << "time mismatch! " << i << ": " 
					<< ts[i] << " != " << firstDay << " or "
					<< ts[i+1] << " != " << lastDay << std::endl;
		}
		for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
			int g = data.weekly_ili[w].gf_by_age[a].first;
			double DeltaS = std::max(0.0, Ss[i][a] - Ss[i+1][a]);
			double bg_se = seasonalBgEffect.getValue();
			double bgili = calcBackgroundIli(bg[a]+bg_se, ampBg, phaseBg, firstDay);
			double lambda = g * q[a] * (DeltaS + bgili - DeltaS*bgili);
			int f = rng.polya(lambda, dispersion); // simulate GP visits
			sdata.weekly_ili[w].gf_by_age[a].second = f;
		}
		// todo: special case for the age aggregated data?
	}
	
	// Consultation
	
	if ( data.consultation.observed ) {
		for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
			int f = data.consultation.fd_by_age[a].first;
			int d = rng.binomial(f, q[a]);
			sdata.consultation.fd_by_age[a].second = d;
		}
	}
	
	return sdata;
}

void McmcSeasonState::update(double T, // temperature for path sampling
		const double bg[], const double q[], // noise and reporting
		double dispersion,
		double beta, // epidemic parameters
		double ampBg, double phaseBg, // oscilation of bg ILI
		const double meanLogitS0[], const double sdLogitS0[], // hyper parameters
		double meanT0, double sdT0, // hyper parameters
		double meanLogGamma, double sdLogGamma, // recovery rate hypers
		double meanLogKappa, double sdLogKappa, // incubation rate hypers
		double sdSeasonalBgEffect, // hyper
		IliDataStruct & data, // data
		Rng & rng) { // random number generator
	// t0...
	updateSinglePar(T, bg, q, dispersion, beta, ampBg, phaseBg, 
			meanLogitS0, sdLogitS0, meanT0, sdT0, 
			meanLogGamma, sdLogGamma, meanLogKappa, sdLogKappa,
			sdSeasonalBgEffect,
			data, t0, rng);
	// gamma
	updateSinglePar(T, bg, q, dispersion, beta, ampBg, phaseBg, 
			meanLogitS0, sdLogitS0, meanT0, sdT0, 
			meanLogGamma, sdLogGamma, meanLogKappa, sdLogKappa,
			sdSeasonalBgEffect,
			data, gamma, rng);
	// kappa
	updateSinglePar(T, bg, q, dispersion, beta, ampBg, phaseBg, 
			meanLogitS0, sdLogitS0, meanT0, sdT0, 
			meanLogGamma, sdLogGamma, meanLogKappa, sdLogKappa,
			sdSeasonalBgEffect,
			data, kappa, rng);
	// seasonalBgEffect...
	updateSinglePar(T, bg, q, dispersion, beta, ampBg, phaseBg, 
			meanLogitS0, sdLogitS0, meanT0, sdT0, 
			meanLogGamma, sdLogGamma, meanLogKappa, sdLogKappa,
			sdSeasonalBgEffect,
			data, seasonalBgEffect, rng);
	// S0...
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		updateSinglePar(T, bg, q, dispersion, beta, ampBg, phaseBg, 
				meanLogitS0, sdLogitS0, meanT0, sdT0, 
				meanLogGamma, sdLogGamma, meanLogKappa, sdLogKappa,
				sdSeasonalBgEffect,
				data, S0[a], rng);
	}
}

void McmcSeasonState::updateSinglePar(double T,
			const double bg[], const double q[],
			double dispersion,
			double beta,
			double ampBg, double phaseBg, // oscilation of bg ILI
			const double meanLogitS0[], const double sdLogitS0[],
			double meanT0, double sdT0,
			double meanLogGamma, double sdLogGamma,
			double meanLogKappa, double sdLogKappa,
			double sdSeasonalBgEffect,
			IliDataStruct & data,
			ParPrior & par, Rng & rng) { // the parameter and a RNG
	if ( !par.isLocked() ) {
		// Ss needs to be up-to-date here!!
		double pll_old = posterior_loglike(T, bg, q, dispersion, beta, ampBg, phaseBg,
				meanLogitS0, sdLogitS0, meanT0, sdT0, 
				meanLogGamma, sdLogGamma, meanLogKappa, sdLogKappa,
				sdSeasonalBgEffect, data, false);
		par.mutate(rng);
		/* re-integrate for all SingleSeasonState parameters. Use the other Ss cache
		 * first find out what Ss cache to use (SsA or SsB)
		 */
		double** Ss_old = Ss;
		double** Ss_new = (Ss_old == SsA ? SsB : SsA);
		// then set Ss to the right cache
		Ss = Ss_new;
		// now re-integrate (using Ss == Ss_new)
		double pll_new = posterior_loglike(T, bg, q, dispersion, beta, ampBg, phaseBg,
				meanLogitS0, sdLogitS0, meanT0, sdT0, 
				meanLogGamma, sdLogGamma, meanLogKappa, sdLogKappa,
				sdSeasonalBgEffect, data, true);
		double logu = log(rng.uniform());
		if ( pll_new >= pll_old + logu ) {
			par.ac1up(); // accept!
			// Ss can stay equal to Ss_new, the other cache (Ss_old) is invalid
		}
		else { 
			par.go_back(); // reject!
			// Ss_new is invalid, set Ss back to Ss_old
			Ss = Ss_old;
			// manually reset ll_cache and pll_cache
			resetLoglikeCache();
			resetPostLoglikeCache();
		}
		if ( mode == WARMUP_MODE ) {
			par.updatePvar(); // only optimize pvar if in WARMUP mode...
		}
	} // if ( !par.isLocked() ) else do nothing
}

void McmcSeasonState::print(std::ostream & os) const {
    os << seasonString(seasonNr) << " ";
	os << ( !ll_cache_good ? BRED_TERM_FONT : "" ) << "ll: " << ll_cache << " " << DEF_TERM_FONT;
	os << ( R0 <= 1.0 ? BRED_TERM_FONT : "" ) << "R_0: " << R0 << " " << DEF_TERM_FONT;
	double lengthOfEpidemic = (tEnd.getValue() - tStart.getValue())/NUMBER_OF_WEEKDAYS;
	os << ( !tEndDefined || !tStartDefined ? BRED_TERM_FONT : "" ) 
	   << "length: " << lengthOfEpidemic << " " << DEF_TERM_FONT;
	os << "t_0: " << t0 << " ";
	os << "gm.: " << gamma << " ";
	os << "kp.: " << kappa << " ";
	os << "bg: " << seasonalBgEffect << " ";
	os << "S_0:";
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		os << " " << S0[a];
	}
}

std::string McmcSeasonState::xmlString() const {
	std::stringstream ss;
	ss << "<mcmc_ss_state "
	   << "ll='" << ll_cache << "' "
  	   << "ll_valid='" << (ll_cache_good ? 1 : 0) << "' "
   	   << "pll='" << pll_cache << "' "
  	   << "pll_valid='" << (pll_cache_good ? 1 : 0) << "' "
	   << "season='" << seasonNr << "' "
	   << "t_start_valid='" << (tStartDefined ? 1 : 0) << "' " 
  	   << "t_end_valid='" << (tEndDefined ? 1 : 0) << "' " 
	   << ">" << std::endl;
	// xml strings of parameters
	ss << t0.xmlString() << std::endl;
	ss << gamma.xmlString() << std::endl;
	ss << kappa.xmlString() << std::endl;
	ss << seasonalBgEffect.xmlString() << std::endl;
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		ss << S0[a].xmlString() << std::endl;
	}
	ss << R0.xmlString() << std::endl;
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		ss << Sinf[a].xmlString() << std::endl;
	}
	ss << tStart.xmlString() << std::endl;
	ss << tEnd.xmlString() << std::endl;
	// other stuff??
	ss << "</mcmc_ss_state>";
	return ss.str();
}

std::ostream & operator<<(std::ostream & os, const McmcSeasonState & state) {
	state.print(os);
	return os;
}

int McmcSeasonState::countParameters() const {
	int k = 0;
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		if ( !S0[a].isLocked() ) ++k;
	}
	if ( !t0.isLocked() ) ++k;
	if ( !gamma.isLocked() ) ++k;
	if ( !kappa.isLocked() ) ++k;
	if ( !seasonalBgEffect.isLocked() ) ++k;
	return k;
}

/** methods for SsUpdateJob, used by McmcFullState::update **/

SsUpdateJob::SsUpdateJob(McmcSeasonState* seasonState, int s, double T,
		double* meanBg, double* slopeBg, 
		double* offset, double* slope, // i.e. reporting (todo: rename)
		double dispersion, double beta, 
		double ampBg, double phaseBg, 
		double* meanLogitS0, double* sdLogitS0, double* slopeLogitS0, 
		double meanT0, double sdT0, 
		double meanLogGamma, double sdLogGamma,
		double meanLogKappa, double sdLogKappa,
		double sdSeasonalBgEffect,
		IliDataStruct* data) :
		seasonState(seasonState), s(s), T(T), 
		meanBg(meanBg), slopeBg(slopeBg), offset(offset), slope(slope),
		dispersion(dispersion), beta(beta), ampBg(ampBg), phaseBg(phaseBg),
		meanLogitS0(meanLogitS0), sdLogitS0(sdLogitS0), slopeLogitS0(slopeLogitS0), 
		meanT0(meanT0), sdT0(sdT0), meanLogGamma(meanLogGamma), sdLogGamma(sdLogGamma),
		meanLogKappa(meanLogKappa), sdLogKappa(sdLogKappa),
		sdSeasonalBgEffect(sdSeasonalBgEffect),
		data(data) {
	/* empty */
}

bool SsUpdateJob::execute(Rng & tlRng) {
	// get values of parameters that are season specific
	double q[NUMBER_OF_AGE_CLASSES];
	double bg[NUMBER_OF_AGE_CLASSES];
	double meanLogitS0vals[NUMBER_OF_AGE_CLASSES];
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		bg[a] = calcLinearTrend(meanBg[a], slopeBg[a], s, CENTRAL_SEASON);
		q[a] = calcTrendendProb(offset[a], slope[a], s, CENTRAL_SEASON);
		meanLogitS0vals[a] = calcLinearTrend(meanLogitS0[a], slopeLogitS0[a], s, CENTRAL_SEASON);
	}
	// the update the ss state
	seasonState->update(T, bg, q, dispersion, beta, ampBg, phaseBg,
			meanLogitS0vals, sdLogitS0, meanT0, sdT0, 
			meanLogGamma, sdLogGamma, meanLogKappa, sdLogKappa,
			sdSeasonalBgEffect, 
			*data, tlRng);
	return true;
}

/** methods for SsLoglikeJob, used by McmcFullState::loglike **/

SsLoglikeJob::SsLoglikeJob(McmcSeasonState* seasonState, int s,
		double* meanBg,	double* slopeBg,
		double* offset, double* slope, // reporting
		double dispersion, double beta,
		double ampBg, double phaseBg,
		bool reintegrate, bool record_lls, 
		IliDataStruct* data) :
		seasonState(seasonState), s(s), meanBg(meanBg), slopeBg(slopeBg),
		offset(offset), slope(slope), dispersion(dispersion), beta(beta),
		ampBg(ampBg), phaseBg(phaseBg), 
		reintegrate(reintegrate), record_lls(record_lls),
		data(data) {
	// empty
}

bool SsLoglikeJob::execute(Rng & tlRng) {
	// get values of parameters that are season specific
	double q[NUMBER_OF_AGE_CLASSES];
	double bg[NUMBER_OF_AGE_CLASSES];
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		bg[a] = calcLinearTrend(meanBg[a], slopeBg[a], s, CENTRAL_SEASON);
		q[a] = calcTrendendProb(offset[a], slope[a], s, CENTRAL_SEASON);
	}
	// the loglike member function caches the computed log likelihood
	seasonState->loglike(bg, q, dispersion, beta, ampBg, phaseBg,
			*data, reintegrate, record_lls);
	return true;
}
