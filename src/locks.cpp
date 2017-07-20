#include "locks.hpp"

void Lock::print(std::ostream & os) const {
	/** Lock::print could be used to print 
	 * the entries of a lock file
	 * TODO: write function xmlString.
	 */
	switch ( type ) {
		case INVALID: {
			os << "# WARNING: invalid lock";
			break;
		}
		case CONSTANT: {
			os << "CONSTANT" << "\t" 
			   << name       << "\t " 
			   << target_value;
			break;
		}
		case EQUALS: {
			os << "EQUALS" << "\t"
			   << name     << "\t"
			   << target_name;
			break;
		}
	}
}

std::ostream & operator<<(std::ostream & os, const Lock & lock) {
	lock.print(os); return os;
}



bool importLocks(std::string filename, std::list<Lock> & locks) {
	/** locks is a (possibly) empty list, passed by ref and will be 
	 * "filled: by this function.
	 * The file must contain tab separated lines with fields indicating
	 * 1) lock type (CONSTANT or EQUALS)
	 * 2) locked parameter name
	 * 3) target of the lock (can be a double or a string)
	 */
	// some aux data for parsing...
	enum LockFields { TYPE_FIELD=0, NAME_FIELD, TARGET_FIELD };
	std::map<std::string, Lock::LockType> typeMap;
	typeMap["CONSTANT"] = Lock::CONSTANT;
	typeMap["EQUALS"] = Lock::EQUALS;
	/* todo: in the future, there may be others... 
	 * Smarter way of defining the typeMap?
	 */
	
	bool succes = true;
	std::fstream fileHandle;
	fileHandle.open(filename.c_str(), std::fstream::in);
	if ( !fileHandle.good() ) succes = false;
	while ( fileHandle.good() ) {
		std::string line;
		std::getline(fileHandle, line);
		// skip comments and empty lines
		if ( line.size() == 0 || line[0] == COMMENT_CHARACTER ) {
			continue;
		}
		std::istringstream lineHandle(line);
		Lock lock;
		int field = TYPE_FIELD;
		while ( lineHandle.good() ) {
			std::string word;
			std::getline(lineHandle, word, '\t');
			switch ( field ) {
				case TYPE_FIELD: {
					lock.type = typeMap[word];
					break;
				}
				case NAME_FIELD: {
					lock.name = word;
					break;
				}
				case TARGET_FIELD: {
					// behaviour depends on lock.type...
					switch ( lock.type ) {
						case Lock::INVALID: {
							// todo: warning?
							break;
						}
						case Lock::CONSTANT: {
							std::stringstream(word) >> lock.target_value;
							break;
						}
						case Lock::EQUALS: {
							lock.target_name = word;
							break;
						}
						default: {
							// todo: warning?
						}
					}
					break;
				}
				default: {
					// todo: warning?
					break;
				}
			}
			++field;
		}
		locks.push_back(lock);
	}
	return succes;
}
