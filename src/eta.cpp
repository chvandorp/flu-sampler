#include "eta.hpp"

EtaEstimator::EtaEstimator(int N) :
        ct(0.0), etl(0.0), n(0), N(N) {
    tic = time(NULL);
    toc = tic;
}
 
void EtaEstimator::update() {
    toc = time(NULL); // the current time
    ct = difftime(toc, tic); // time between start and now
    ++n; // number of steps taken
    etl = std::max((ct/n) * (N-n), 0.0); // expected time left
}
 
void EtaEstimator::print(std::ostream & os) const {
    double etlprime = etl;
    int days = floor(etlprime / secperday);
    etlprime -= days * secperday;
    int hours = floor(etlprime / secperhour); 
    etlprime -= hours * secperhour;
    int minutes = floor(etlprime / secperminute);
    etlprime -= minutes * secperminute;
    int seconds = floor(etlprime);
    os << (days > 0 ? std::to_string(days) + " " : "")
       << hours << ":"
       << (minutes < 10 ? "0" : "") << minutes << ":"
       << (seconds < 10 ? "0" : "") << seconds;
}
 
std::ostream & operator<<(std::ostream & os,
        const EtaEstimator & eta) {
    eta.print(os);
    return os;
}
