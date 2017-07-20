#include "exceptions.hpp"

const char* MsgException::what() const throw() {
    return ("Exception: " + msg).c_str();
}
