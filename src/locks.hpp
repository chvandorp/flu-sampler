/** locks.cpp/hpp defines a couple of methods concerning the "locking"
 * of parameters in order to select a model.
 * the class ParPrior has a method setLock that prevents the sampler
 * from updating the parameter, effectively making it a contant.
 * ParPrior::setLock can also force a parameter to be equal to another,
 * allowing for (somewhat resticted) relations.
 * locks.cpp/hpp defines the struct Lock and a function to read locks 
 * from a file.
 * The lock file accepts wildcards, that are expanded by
 * McmcState::setLocks(...). The current implementation uses C++ <regex>
 * and e.g. the wildcard * must be written as .*
 * TODO: Maybe we could convert more basic wildcards into C++ regex.
 */

#ifndef LOCKS_HPP_
#define LOCKS_HPP_

#include <string>
#include <list>
#include <map> // switch with strings
#include <iostream>
#include <sstream>
#include <fstream>

#include "constants.hpp" // defines COMMENT_CHARACTER

struct Lock {
	enum LockType { INVALID=0, CONSTANT, EQUALS } type;
	Lock() : type(INVALID) {};
	std::string name;
	std::string target_name;
	double target_value;
	void print(std::ostream & ) const;
};

std::ostream & operator<<(std::ostream & , const Lock & );

bool importLocks(std::string , std::list<Lock> & );

#endif
