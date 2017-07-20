#ifndef EXCEPTIONS_HPP
#define EXCEPTIONS_HPP

#include <exception>
#include <string>

#define RIGHT_HERE (std::string(" (in function ") + __FUNCTION__ + ") ")

class MsgException : public std::exception {
public:
	MsgException(std::string msg="") : msg(msg) { /* empty*/ }
	~MsgException() throw() { /* empty*/ }
	virtual const char* what() const throw();
protected:
	std::string msg;
};

#endif
