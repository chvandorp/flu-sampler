#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <cmath> // fmod, sqrt...
#include <limits> // for infinity

#include "distributions.hpp" // includes Rng
#include "constants.hpp"

/** here we define objects that serve as parameters of the models.
 * In the current algorithm, the parameter class has a pointer 
 * to a prior distribution, and handles the proposal distribution.
 * TODO: 1) make one class for discrete and continuous variables\
 * 2) allow for parameter vectors and multivariate proposals.
 */

enum HomotopyClass {
	CONTRACTIBLE=0,
	CIRCULAR
};

class ParPrior {
public:
	ParPrior();
	ParPrior(double , double , Prior* , std::string name=""); 
	// initial value, proposal variance, a prior, and an optional name
	ParPrior(const ParPrior & ); 
	// copy of prior, uses Prior::dup method
	ParPrior & operator=(const ParPrior & ); 
	// copy of prior, uses Prior::dup method
	~ParPrior(); // deletes prior
	virtual ParPrior* dup() const; // make a copy and return pointer
	void lock(double ); // keep constant to argument
	void lock(ParPrior* ); // keep equal to argument
	bool isLocked() const;
	void setBounds(double , double );
	void setLBound(double );
	void setUBound(double );
    bool isBounded() const; // uboundbool && lboundbool
    double getLengthInterval() const; // inf if unbounded...
	void setHomotopy(HomotopyClass );
	double getValue() const;
	double getPvar() const;
	double getAr() const; // total acceptance ratio
	double loglike(Transformation* fun=NULL) const;
	void mutate(Rng & );
	void ac1up();
	void go_back();
	void updatePvar();
	void resetUpdateCounter(); // sets uc to 0: allows pvar to be re-optimized.
	std::string getName() const;
	std::string xmlString() const;
	void print(std::ostream & ) const;
protected:
	std::string name;
	double value;
	double old_value;
	double pvar;
	bool lboundbool;
	double lbound;
	bool uboundbool;
	double ubound;
	HomotopyClass homotopy_support; // contractible, circular, ...
	Prior* prior;
	int stepcounter; // counts number of steps taken
	int acceptcounter; // overall
	int ac; // acceptance counter (short time scale)
	int uc; // update counter
	bool locked; // if true, then the parameter is constant or equal to another...
	ParPrior* lockedTo; // if not NULL, the value is equal to lockedTo->getValue()
};

std::ostream & operator<<(std::ostream & , const ParPrior & );


class CompoundPar {
public:
	CompoundPar() : value(0.0) {}
	CompoundPar(double value, std::string name="") : name(name), value(value) {}
	CompoundPar & operator=(double );
	bool operator<(double ) const;
	bool operator>(double ) const;
	bool operator==(double ) const;
	bool operator<=(double ) const;
	bool operator>=(double ) const;
	void setValue(double );
	double getValue() const;
	std::string getName() const;
	std::string xmlString() const;
	void print(std::ostream & ) const;
protected:
	std::string name;
	double value;
};

std::ostream & operator<<(std::ostream & , const CompoundPar & );


/** Simple class for discrete parameters. UNDER CONSTRUCRTION!
 */

class DiscreteParPrior {
public:
	DiscreteParPrior();
	DiscreteParPrior(int , double , std::string name="");
	void setBounds(int , int );
	void setLBound(int );
	void setUBound(int );
	void lock(int ); // keep constant to argument
	bool isLocked() const;
	int getValue() const;
	void mutate(Rng & );
	void go_back();
protected:
	std::string name;
	int value;
	int old_value;
	double jump_prob;
	int lbound, ubound;
	bool lboundbool, uboundbool;
	bool locked;
};

#endif
