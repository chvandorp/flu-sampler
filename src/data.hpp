#ifndef DATA_HPP
#define DATA_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <list>

#include "constants.hpp"
#include "aux.hpp" // for waic export

/* data used in Mcmc */

struct IliWeekDataStruct {
	IliWeekDataStruct();
	int firstDay;
	int lastDay;
	std::pair<int, int> gf_by_age[NUMBER_OF_AGE_CLASSES];
	std::pair<int, int> gf_total;
	bool age_stratified;
	std::string xmlString() const;
	std::string singleLinePrint() const;
	void print(std::ostream & ) const;
	int countObservations() const;
	// members used for WAIC computation	
	std::list<double> lls_by_age[NUMBER_OF_AGE_CLASSES];
	std::list<double> lls_total;
	void clear_lls(); // sometimes we dont want to remember all lls
	std::string waicXmlString() const;
};

std::ostream & operator<<(std::ostream & , const IliWeekDataStruct & );

bool compareIliWeekData(const IliWeekDataStruct & , 
	const IliWeekDataStruct & );

struct ConsultDataStruct {
	ConsultDataStruct();
	int seasonNr;
    // yearly info on the number of ILI+ indidividuals consulting a GP
    std::pair<int, int> fd_by_age[NUMBER_OF_AGE_CLASSES];
    std::pair<int, int> fd_total; // todo
    bool observed; // consultation data is not available for many seasons
    bool age_stratified; // todo
    std::string xmlString() const;
   	void print(std::ostream & ) const;
   	int countObservations() const;
    std::list<double> lls_by_age[NUMBER_OF_AGE_CLASSES];
    std::list<double> lls_total;
   	void clear_lls(); // sometimes we dont want to remember all lls
    std::string waicXmlString() const;
};

std::ostream & operator<<(std::ostream & , const ConsultDataStruct & );

struct IliDataStruct {
	IliDataStruct();
	~IliDataStruct(); // delete[] ts and C
	int seasonNr;
	int numberOfWeeks;
	int firstDay;
	int lastDay;
	std::vector<IliWeekDataStruct> weekly_ili;
	ConsultDataStruct consultation;
	double* ts; // timepoints: useful to remember for integrator
	int numberUniqueTimePoints; // length of ts
	double* C; // this seasons contact matrix
    // methods for printing
	std::string xmlString() const;
	std::string singleLinePrint() const;
	void print(std::ostream & ) const;
	int countObservations() const;
	void clear_lls();
	std::string waicXmlString() const;
};

std::ostream & operator<<(std::ostream & , const IliDataStruct & );
std::ostream & operator<<(std::ostream & , const std::vector<IliDataStruct> & );

/* data import functions... These functions are complicated,
 * and all do rather similar things. Is there a better way?
 */

bool importContactMatrix(std::string, double* );
bool importIliData(std::string, std::vector<IliDataStruct> & );
bool importTimingSeasons(std::string, std::vector<IliDataStruct> & );
bool importGisData(std::string, std::vector<IliDataStruct> & );

bool addContactMatrixToIliData(double*, std::vector<IliDataStruct> & ); 

bool readIliLine(std::string line , IliWeekDataStruct & , int & );

int countIliObservations(const std::vector<IliDataStruct> & ); // for BIC and WBIC

#endif
