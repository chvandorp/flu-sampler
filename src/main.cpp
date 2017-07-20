#include "main.hpp"

int main(int argc, char** argv) {
	std::cout << "# starting flu-sampler (version 2.1)..." << std::endl;
 
	std::string mode = "help";
    std::string vid = "test";
    unsigned long seed = 144169;
    int length = 1000;
    int thinning = 10;
    int imitate = 20; // need to get 2.5 and 97.5 percentiles

	// ignore all other input when the help flag (-h) is given
	bool foundHelpFlag = false;
	int c = 0;
    while ( (c = getopt(argc, argv, "m:v:s:l:i:t:h")) != -1 && !foundHelpFlag ) {
        switch (c) {
			case 'm': { // mode
				std::stringstream ss;
				ss << optarg;
				mode = ss.str();
				break;
			}
			case 'v': { // version id
				std::stringstream ss;
				ss << optarg;
				vid = ss.str();
				break;
			}
			case 's': { // seed
				std::stringstream ss;
				ss << optarg;
				ss >> seed;
				break;
			}
			case 'l': { // length of chain
				std::stringstream ss;
				ss << optarg;
				ss >> length;
				break;
			}
			case 't': { // thinning
				std::stringstream ss;
				ss << optarg;
				ss >> thinning;
				break;
			}
			case 'i': { // number of simulations
				std::stringstream ss;
				ss << optarg;
				ss >> imitate;
                break;
			}
			case 'h': { // force help mode
				foundHelpFlag = true;
				mode = "help";
				break;
			}
			case '?': {
				break;
			}
			default: {
				std::stringstream ss; 
				ss << "getopt returned nonsense" << RIGHT_HERE;
				throw MsgException(ss.str());
			}
        }
    }

	// define modes
    enum ModeId { HELP_MODE=1, MCMC_MODE, WBIC_MODE, BFAC_MODE }; // skip 0...
    std::map<std::string, int> modeMap;
    modeMap["help"] = HELP_MODE;
    modeMap["mcmc"] = MCMC_MODE;
    modeMap["wbic"] = WBIC_MODE;
            
    // print input for verification
	std::cout << "# version id: " << vid  << std::endl
	          << "# seed:       " << seed << std::endl
	          << "# length:     " << length << std::endl
	          << "# thinning:   " << thinning << std::endl
	          << "# imitations: " << imitate << std::endl
	          << "# mode:       " << mode << std::endl;
	          
	// give the user the chance to ctrl-c if they made a mistake
	if ( modeMap[mode] != HELP_MODE ) {
		std::cout << "> ok? (press enter)";
		std::string str;
		std::getline(std::cin, str); // accepts characters until return
	}
		
    // select mode    
    switch ( modeMap[mode] ) { // assumes that modeMap["unknown string"] is always 0
		case HELP_MODE:
			printHelp();
			break;
		case MCMC_MODE:
			runMcmc(vid, seed, length, thinning, imitate, false);
			break;
		case WBIC_MODE:
			runMcmc(vid, seed, length, thinning, imitate, true);
			break;
		default:
			std::cout << "/!\\ invalid mode given..." << std::endl;
			printHelp();
			break;
	}		
	std::cout << "# goodbye!" << std::endl;
	return 0;
}

void runMcmc(std::string vid, unsigned long seed, int length, 
             int thinning, int imitate, bool wbic_mode) {
	// parameters (todo: get from command line)
	int burnin = length; // cf. Stan default
	
	std::cout << "# starting MCMC estimation..." << std::endl;
	// make an instance of the Mcmc object
	Mcmc chain(length, burnin, thinning, wbic_mode, seed, vid);
	chain.run();
	
	std::cout << "# writing chain to file..." << std::endl;
	// write to xml file (TODO: put this printing somewhere else)
	std::string chainXmlFileName = "data/full-ili-chain-" + vid + ".xml";
	std::fstream chainXmlFileHandle;
	chainXmlFileHandle.open(chainXmlFileName.c_str(), std::fstream::out);
	chainXmlFileHandle << chain;
	chainXmlFileHandle.close();
	
	std::cout << "# writing pointwise waic to file..." << std::endl;
	std::string waicXmlFileName = "data/pointwise-waic-" + vid + ".xml";
	std::fstream waicXmlFileHandle;
	waicXmlFileHandle.open(waicXmlFileName.c_str(), std::fstream::out);
	chain.printPointwiseWaic(waicXmlFileHandle);
	waicXmlFileHandle.close();

	std::cout << "# sampling from posterior and simulating data..." << std::endl;
	// sample and simulate
	std::string simXmlFileName = "data/sili-" + vid + ".xml";
	std::fstream simXmlFileHandle;
	simXmlFileHandle.open(simXmlFileName.c_str(), std::fstream::out);
	chain.simulateAndPrint(simXmlFileHandle, imitate);
	simXmlFileHandle.close();
	
	// chain.runSingleSeason(43); // testing...
	
	// TODO: throw exceptions when the files cannot be opened
}

void printHelp() {
	/** a help message is stored in doc/helpmsg.txt
	 * Print this help message and return
	 */
	std::fstream fileHandle;
	fileHandle.open("doc/helpmsg.txt", std::fstream::in);
	if ( fileHandle.good() ) {
		while ( fileHandle.good() ) {
			std::string line;
			std::getline(fileHandle, line);
			std::cout << line << std::endl;
		}
		fileHandle.close();
	}
	else {
		std::cerr << "ERROR: can't find help message... "
		             "You're on your own :-(" << std::endl;
	}
}
