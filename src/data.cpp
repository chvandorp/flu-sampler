#include "data.hpp"

/** methods for structs: */

/** Methods for IliWeekDataStruct */

IliWeekDataStruct::IliWeekDataStruct() {
	firstDay = 0;
	lastDay = NUMBER_OF_WEEKDAYS;
	age_stratified = false;
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		gf_by_age[a] = std::pair<int, int>(0, 0);
	}
	gf_total = std::pair<int, int>(0, 0);
}

std::string IliWeekDataStruct::singleLinePrint() const {
	std::stringstream ss;
	ss << firstDay;
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		ss << " " << gf_by_age[a].first << " " << gf_by_age[a].second;
	}
	return ss.str();
}

std::string IliWeekDataStruct::xmlString() const {
	std::stringstream ss;
	ss << "<ili_week_data "
	   << "first_day='" << firstDay << "' "
	   << "agestrat='" << ( age_stratified ? "true" : "false" ) << "' "
	   << ">" << std::endl;
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		ss << "<ili_gf_data "
		   << "ageclass='" << a << "' "
		   << "g='" << gf_by_age[a].first << "' " // GP population
		   << "f='" << gf_by_age[a].second << "' " // with flu-like illness
		   << "/>" << std::endl;
	}
	ss << "</ili_week_data>"; // no endl!
	// TODO: add gf_total (and more?)
	return ss.str();
}

void IliWeekDataStruct::print(std::ostream & os) const {
	// os << singleLinePrint();
	os << xmlString();
}

int IliWeekDataStruct::countObservations() const {
	if ( age_stratified ) {
		return NUMBER_OF_AGE_CLASSES;
	}
	else { // only one observation (fraction Xia)
		return 1;
	}
}

void IliWeekDataStruct::clear_lls() {
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		lls_by_age[a].clear();
	}
	lls_total.clear();
}

std::string IliWeekDataStruct::waicXmlString() const {
	std::stringstream ss;
	ss << "<waic_week_data "
	   << "first_day='" << firstDay << "' "
	   << "agestrat='" << ( age_stratified ? "true" : "false" ) << "' "
	   << ">" << std::endl;
	if ( age_stratified ) {
		for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
			ss << "<waic_gf_data "
			   << "ageclass='" << a << "' "
			   << "lpd='" << log_mean_exp(lls_by_age[a]) << "' "
			   << "pwaic='" << sample_variance(lls_by_age[a]) << "' "
			   << "/>" << std::endl;
		}
	}
	else {
		ss << "<waic_gf_data "
		   << "lpd='" << log_mean_exp(lls_total) << "' "
		   << "pwaic='" << sample_variance(lls_total) << "' "
		   << "/>" << std::endl;		
	}
	ss << "</waic_week_data>";
	return ss.str();
}

std::ostream & operator<<(std::ostream & os, const IliWeekDataStruct & weekly_ili) {
	weekly_ili.print(os);
	return os;
}


/** Methods for ConsultDataStruct */

ConsultDataStruct::ConsultDataStruct() {
	observed = false;
	age_stratified = false;
	seasonNr = 0;
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		fd_by_age[a] = std::pair<int, int>(0, 0);
	}
	fd_total = std::pair<int, int>(0, 0);
}

std::string ConsultDataStruct::xmlString() const {
	std::stringstream ss;
	ss << "<consult_data "
	   << "season='" << seasonNr << "' "
	   << "observed='" << ( observed ? "true" : "false" ) << "' "
   	   << "agestrat='" << ( age_stratified ? "true" : "false" ) << "' "
	   << ">" << std::endl;
	if ( observed ) {
		for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
			ss << "<consult_fd_data "
			   << "ageclass='" << a << "' "
		       << "f='" << fd_by_age[a].first << "' " // flu-like illness
		       << "d='" << fd_by_age[a].second << "' " // GP consult
		       << "/>" << std::endl;
		}
	}
	ss << "</consult_data>"; // no endl!
	// TODO: add fd_total (and more?)
	return ss.str();
}

int ConsultDataStruct::countObservations() const {
	if ( observed ) {
		if ( age_stratified ) {
			return NUMBER_OF_AGE_CLASSES;
		}
		else { // at the moment redundant
			return 1;
		}
	} // zero observations...
	else {
		return 0;
	}
}

void ConsultDataStruct::clear_lls() {
	for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
		lls_by_age[a].clear();
	}
	lls_total.clear();
}

std::string ConsultDataStruct::waicXmlString() const {
	std::stringstream ss;
	ss << "<waic_consult_data "
	   << "season='" << seasonNr << "' "
	   << "observed='" << ( observed ? "true" : "false" ) << "' "
   	   << "agestrat='" << ( age_stratified ? "true" : "false" ) << "' "
	   << ">" << std::endl;
	if ( observed ) {
		if ( age_stratified ) {
			for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
				ss << "<waic_consult_fd_data "
			       << "ageclass='" << a << "' "
			       << "lpd='" << log_mean_exp(lls_by_age[a]) << "' "
			       << "pwaic='" << sample_variance(lls_by_age[a]) << "' "
			       << "/>" << std::endl;
			}
		}
		else {
			ss << "<waic_consult_fd_data "
			   << "lpd='" << log_mean_exp(lls_total) << "' "
			   << "pwaic='" << sample_variance(lls_total) << "' "
			   << "/>" << std::endl;		
		}
	} // if observed
	ss << "</waic_consult_data>";
	return ss.str();
}

void ConsultDataStruct::print(std::ostream & os) const {
	os << xmlString();
}

std::ostream & operator<<(std::ostream & os, const ConsultDataStruct & consultation) {
	consultation.print(os);
	return os;
}

/** methods for IliDataStruct */

IliDataStruct::IliDataStruct() {
	seasonNr = 0;
	numberOfWeeks = 0;
	firstDay = 0;
	lastDay = 0;
	numberUniqueTimePoints = 0;
	ts = NULL;
	C = NULL;
}

IliDataStruct::~IliDataStruct() {
	delete[] ts;
	delete[] C;
}

std::string IliDataStruct::xmlString() const {
	std::stringstream ss;
	ss << "<ili_season_data "
	   << "season='" << seasonNr << "' "
	   << "first_day='" << firstDay << "' "
	   << "last_day='" << lastDay << "' "
	   << "num_time_points='" << numberUniqueTimePoints << "' " 
	   // for redundancy (TODO why not numberOfWeeks??)
	   << ">" << std::endl;
	for ( int w = 0; w < numberOfWeeks; ++w ) {
		ss << weekly_ili[w] << std::endl;
	}
	ss << consultation << std::endl;
	ss << "</ili_season_data>"; // no endl!
	return ss.str();
}

std::string IliDataStruct::singleLinePrint() const {
	std::stringstream ss;
	for ( int w = 1; w <= numberOfWeeks; ++w ) {
		ss << weekly_ili[w].singleLinePrint() << (w==numberOfWeeks ? "" : " ");
	}
	return ss.str();
}

void IliDataStruct::print(std::ostream & os) const {
	// os << singleLinePrint();
	os << xmlString();
}

int IliDataStruct::countObservations() const {
	int n = 0;
	for ( int w = 0; w < numberOfWeeks; ++w ) {
		n += weekly_ili[w].countObservations();
	}
	n += consultation.countObservations();
	return n;
}

void IliDataStruct::clear_lls() {
    // clear lls for weekly ILI data
	for ( int w = 0; w < numberOfWeeks; ++w ) {
		weekly_ili[w].clear_lls();
	}    
    // clear seasonal lls
    for ( int a = 0; a < NUMBER_OF_AGE_CLASSES; ++a ) {
        consultation.lls_by_age[a].clear();
    }
    consultation.lls_total.clear();
}

std::string IliDataStruct::waicXmlString() const {
	std::stringstream ss;
	ss << "<waic_season_data "
       << "season='" << seasonNr << "' "
       << ">" << std::endl;
	for ( int w = 0; w < numberOfWeeks; ++w ) {
		ss << weekly_ili[w].waicXmlString() << std::endl;
	}
	ss << consultation.xmlString() << std::endl;
	ss << "</waic_season_data>";
	return ss.str();
}


std::ostream & operator<<(std::ostream & os, const IliDataStruct & seasonal_ili) {
	seasonal_ili.print(os);
	return os;
}

std::ostream & operator<<(std::ostream & os, const std::vector<IliDataStruct> & ili) {
	os << "<ili_data >" << std::endl; // todo: add parameters
	std::vector<IliDataStruct>::const_iterator vit = ili.begin();
	for ( ; vit != ili.end(); ++vit ) {
		os << (*vit) << std::endl;
	}
	os << "</ili_data>";
	return os;
}

bool compareIliWeekData(const IliWeekDataStruct & weekly_ili1, 
		const IliWeekDataStruct & weekly_ili2) {
	if ( weekly_ili1.firstDay < weekly_ili2.firstDay ) return true;
	else {
		if ( weekly_ili1.firstDay == weekly_ili2.firstDay ) {
			return weekly_ili1.lastDay < weekly_ili2.lastDay;
		}
		else return false;
	}
}



/* functions for importing data: */

bool importContactMatrix(std::string fileName, double* C) {
	bool status = true;
	
	std::fstream fileHandle;
	fileHandle.open(fileName.c_str(), std::fstream::in);

	if ( fileHandle.good() ) {
		std::string line;
		int i = 0;
		while ( fileHandle.good() && i < NUMBER_OF_AGE_CLASSES ) {
			std::getline(fileHandle, line);
			std::stringstream lineHandle(line);
			std::string word;
			int j = 0;
			while ( lineHandle.good() && j < NUMBER_OF_AGE_CLASSES ) {
				std::getline(lineHandle, word, '\t');
				double val = 0.0;
				std::stringstream(word) >> val;
				C[NUMBER_OF_AGE_CLASSES*i + j] = val; 
				// fill matrix by row (gsl_matrix_view_array interpretation)
				++j;
			}
			++i;
		}
		fileHandle.close();
	}
	else { 
		status = false;
	}
	return status;
}

bool importIliData(std::string fileName, std::vector<IliDataStruct> & data) {
	// create the data structure
	data.resize(NUMBER_OF_SEASONS);
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		data[s].seasonNr = s;
	}
	// now read actual data...
	bool status = true;
	
	std::fstream fileHandle;
	fileHandle.open(fileName.c_str(), std::fstream::in);

	if ( fileHandle.good() ) {
		// first scan through data to get some "meta-data"
		std::string line;
		int linenumber = 1; // for error handling. start counting with 1 (text editors do)
		while ( fileHandle.good() ) {
			std::getline(fileHandle, line);
			IliWeekDataStruct weekly_ili;
			int season;
			// weekly_ili and season passed by ref!
			if ( line.size() > 0 && line[0] != COMMENT_CHARACTER ) {
				if ( readIliLine(line, weekly_ili, season) ) {
					data[season].numberOfWeeks++;
				}
				else {
					std::cerr << "faulty ILI line (" << linenumber << "): " << line << std::endl;
				}
			}
			linenumber++;
		}
		// modify the data structure
		for ( int season = 0; season < NUMBER_OF_SEASONS; ++season ) {
			int numberOfWeeks = data[season].numberOfWeeks;
			data[season].weekly_ili.reserve(numberOfWeeks);
		}
		// now load the real data
		fileHandle.clear(); // re-set eof flag
		fileHandle.seekg(0, std::fstream::beg); // jump to beginning
		while ( fileHandle.good() ) {
			std::getline(fileHandle, line);
			// weekly_ili and season passed by ref!
			if ( line.size() > 0 && line[0] != COMMENT_CHARACTER ) {
				IliWeekDataStruct weekly_ili;
				int season;
				if ( readIliLine(line, weekly_ili, season) ) {
					data[season].weekly_ili.push_back(weekly_ili);
				}
			}
		}
		fileHandle.close();
		// post-processing
		for ( int season = 0; season < NUMBER_OF_SEASONS; ++season ) {
			// sort weeks by first day
			std::sort(data[season].weekly_ili.begin(), 
				data[season].weekly_ili.end(), 
				compareIliWeekData);
			// make a list of all time points needed for the integrator
			std::vector<double> ts_vec;
			ts_vec.reserve(2 * data[season].numberOfWeeks);
			for ( int w = 0; w < data[season].numberOfWeeks; ++w ) {
				ts_vec.push_back(data[season].weekly_ili[w].firstDay);
				ts_vec.push_back(data[season].weekly_ili[w].lastDay);
			}
			std::sort(ts_vec.begin(), ts_vec.end()); // sort (not really needed...)
			std::vector<double>::iterator it = std::unique(ts_vec.begin(), ts_vec.end()); 
			// and remove duplicate points.
			ts_vec.resize(std::distance(ts_vec.begin(), it)); 
			// remove meaningless tail of uniqued ts_vec
			data[season].ts = new double[ts_vec.size()];
			for ( int w = 0; w < int(ts_vec.size()); ++w ) {
				data[season].ts[w] = ts_vec[w];
			}
			data[season].numberUniqueTimePoints = ts_vec.size();
		}
	}
	else {
		status = false; // cannot open file
	}
	// return true on success
	return status;
}

bool importTimingSeasons(std::string fileName, std::vector<IliDataStruct> & data) {
	// assumes that data already has been resized
	std::fstream fileHandle;
	fileHandle.open(fileName.c_str() , std::fstream::in);
	if ( fileHandle.good() ) {
		while ( fileHandle.good() ) {
			std::string line;
			std::getline(fileHandle, line);
			int season;
			int firstDay;
			int lastDay;
			std::stringstream lineHandle(line);
			std::string word;
			int j = 0;
			while ( lineHandle.good() ) {
				std::getline(lineHandle, word, '\t');
				switch ( j ) {
					case 0: { // todo: use an enum
						std::stringstream(word) >> season;
						break;
					}
					case 1: {
						std::stringstream(word) >> firstDay;
						break;
					}
					case 2: {
						std::stringstream(word) >> lastDay;
						break;
					}
					default: {
						break; 
					}
				}
				++j;
			}
			if ( 0 <= season && season < NUMBER_OF_SEASONS ) {
				data[season].firstDay = firstDay;
				data[season].lastDay = lastDay;
			}
		}
		return true;
	}
	else return false;
}


bool addContactMatrixToIliData(double* C, std::vector<IliDataStruct> & data) {
	/* assumes that data has been properly resized
	 * and that C has been allocated
	 */
	/** every season has it's own copy of the contact matrix C.
	 * In a future version, the contact matrix could be dependent on the
	 * age structure in the population
	 */
	bool ok = true;
	for ( auto it = data.begin(); it != data.end(); ++it ) {
		if ( it->C == NULL ) {
			int len = NUMBER_OF_AGE_CLASSES * NUMBER_OF_AGE_CLASSES;
			it->C = new double[len];
			for ( int i = 0; i < len; ++i ) it->C[i] = C[i];
		}
		else {
			ok = false;
		}
	}
	return ok;
}

enum GisDataColumn {
    GIS_YEAR_COL = 0,
    GIS_STRAT_COL,
    GIS_NUMBER_OF_META_COLS
};

// TODO: rename: ConsultData/ConsultationData (don't mention GIS)
bool importGisData(std::string fileName, std::vector<IliDataStruct> & data) {
    bool ok = true;

	std::fstream fileHandle;
	fileHandle.open(fileName.c_str(), std::fstream::in);

	if ( fileHandle.good() ) {
		std::string line;
		int linenumber = 1; // for error handling. start counting with 1 (text editors do)
		while ( fileHandle.good() ) {
			std::getline(fileHandle, line);
			if ( line.size() > 0 && line[0] != COMMENT_CHARACTER ) {
				ConsultDataStruct consultation;
       			int year; // GIS table is by year
       			int season = -1; // map year to corresponding season. Use -1 default
				bool age_stratified = false;
                std::stringstream lineHandle(line);
				std::string word;
				int j = 0; // word count
				while ( lineHandle.good() ) {
					std::getline(lineHandle, word, '\t');
					switch ( j ) { // the j-th word
						case GIS_YEAR_COL: { // the year/season
							std::stringstream(word) >> year;
							// find a corresponding season
							season = yearToSeason(year);
							consultation.seasonNr = season;
							consultation.observed = true;
							break;
						}
						case GIS_STRAT_COL: {
							std::stringstream(word) >> std::boolalpha >> age_stratified;
							consultation.age_stratified = age_stratified;
							break;
						}
						default: {
							int a = (j-GIS_NUMBER_OF_META_COLS)/2;
							if ( a < NUMBER_OF_AGE_CLASSES ) {
								int val = 0;
								std::stringstream(word) >> val;
								int parity = (j-GIS_NUMBER_OF_META_COLS) % 2;
								if ( parity == 0 ) { // f d f d -> f
									consultation.fd_by_age[a].first = val;
									consultation.fd_total.first += val;
								}	
								else { // f d f d -> d
									consultation.fd_by_age[a].second = val;
									consultation.fd_total.second += val;
								}
							}
							else {
								// todo: there might be more data, but for now, ignore it
							}
							break; // redundant (default case)	
						} // default
					} // switch
					++j; // word count
				} // while ( lineHandle.good() )
				if ( j >= GIS_NUMBER_OF_META_COLS + 2*NUMBER_OF_AGE_CLASSES ) {
					ok = ok && true;
				}
				// put the ConsultDataStruct in the right IliDataStruct
				if ( season >= 0 && season < (int) data.size() ) {
					data[season].consultation = consultation;
				}
				else {
					std::cerr << "# WARNING: season " << season << " not in data range. " 
							  << "line " << linenumber << std::endl;
					ok = false;
				}
			} // if ( line.size() > 0 ... )
			linenumber++;
		} // while ( fileHandle.good() )
		fileHandle.close();
    }
    return ok;
}


enum IliDataColumn {
	SEASON_COL = 0,
	DAY_COL,
	DATE_COL,
	STRAT_COL,
	NUMBER_OF_META_COLS // should be last!
};

bool readIliLine(std::string line, IliWeekDataStruct & weekly_ili, int & season) {
	weekly_ili.gf_total.first = 0; weekly_ili.gf_total.second = 0;	
	std::stringstream lineHandle(line);
	std::string word;
	int firstDay;
	bool age_stratified;
	int j = 0;
	while ( lineHandle.good() ) {
		std::getline(lineHandle, word, '\t');
		switch ( j ) { // the j-th word
			case SEASON_COL: { // the season
				std::stringstream(word) >> season;
				break;
			}
			case DAY_COL: { // the day number
				std::stringstream(word) >> firstDay;
				weekly_ili.firstDay = firstDay;
				weekly_ili.lastDay = firstDay + NUMBER_OF_WEEKDAYS;
				break;
			}
			case DATE_COL: {
				// todo...
				break;
			}
			case STRAT_COL: {
				std::stringstream(word) >> std::boolalpha >> age_stratified;
				weekly_ili.age_stratified = age_stratified;
				break;
			}
			default: {
				int a = (j-NUMBER_OF_META_COLS)/2;
				if ( a < NUMBER_OF_AGE_CLASSES ) {
					int val = 0;
					std::stringstream(word) >> val;
					int parity = (j-NUMBER_OF_META_COLS) % 2;
					if ( parity == 0 ) { // g f g f...
						weekly_ili.gf_by_age[a].first = val;
						weekly_ili.gf_total.first += val;
					}
					else {
						weekly_ili.gf_by_age[a].second = val;
						weekly_ili.gf_total.second += val;					
					}
				}
				else {
					// todo: there might be more data, but for now, ignore it
				}
				break;
			}
		}
		++j;
	}
	return j >= 2*NUMBER_OF_AGE_CLASSES + NUMBER_OF_META_COLS;
	// todo: better error checking
}

int countIliObservations(const std::vector<IliDataStruct> & data) { // for (W)BIC
	int n = 0;
	for ( int s = 0; s < NUMBER_OF_SEASONS; ++s ) {
		n += data[s].countObservations();
	}
	return n;
}
