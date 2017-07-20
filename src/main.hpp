#ifndef MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include <unistd.h> // for getopt

#include "mcmc.hpp"

int main(int argc, char** argv);

void runMcmc(std::string , unsigned long , int , int , int , bool );
void printHelp();

#endif
