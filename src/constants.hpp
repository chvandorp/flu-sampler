#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <gsl/gsl_math.h> // math constants
#include <string>

#define NUMBER_OF_SEASONS 45
#define CENTRAL_SEASON (0.5*(NUMBER_OF_SEASONS-1))
#define NUMBER_OF_WEEKDAYS 7
#define DAYS_IN_YEAR 365.24
#define COMMENT_CHARACTER '#' // add some comments to the data files

#define INOCULUM_SIZE 1e-6 // parameter used for seeding the epidemc
#define DURATION_INFECTIOUS_PERIOD 3.0 // days
#define DURATION_EXPOSED_PERIOD 1.0 // days
#define NUMBER_OF_INFECTIOUS_STAGES 2
#define NUMBER_OF_EXPOSED_STAGES 0
#define EPIDEMIC_THRESHOLD 5.1e-4 // definition by NIVEL

#define NUMBER_OF_AGE_CLASSES 6 // set equal to 1 for the aggregated model

#define PVAR_UPDATE_INTERVAL 50
#define OPTIMAL_ACCEPTANCE_RATE 0.44 // 0.44 for Metropolis within Gibbs with normal proposal
#define LOW_ACCEPTANCE_RATE 0.1 // warn when below this value
#define HIGH_ACCEPTANCE_RATE 0.9 // warn when above this value

#define YEAR_ZERO 1970

#define DATA_FOLDER std::string("data/essential/")
#define NUMBER_OF_THREADS 15

// colors
#define BBLUE_TERM_FONT std::string("\033[1;34m")
#define BRED_TERM_FONT std::string("\033[1;31m")
#define BLUE_TERM_FONT std::string("\033[0;34m")
#define RED_TERM_FONT std::string("\033[0;31m")
#define CYAN_TERM_FONT std::string("\033[0;36m")
#define YELLOW_TERM_FONT std::string("\033[0;33m")
#define DEF_TERM_FONT std::string("\033[0m")

#endif
