#ifndef ETA_HPP_
#define ETA_HPP_

#include <ctime>
#include <cmath>
#include <algorithm> // std::max
#include <iostream>

/* first n terms of the geometric series:
 * 1 + alpha + alpha^2 + ... + alpha^{n-1}
 */
inline double geometric_sum(double alpha, int n) {
    double an = pow(alpha, n);
    return (an-1) / (alpha-1);
}

class EtaEstimator {
public:
    EtaEstimator(int N); 
    // constuction starts the clock. Pass the number of steps
    void update();
    void print(std::ostream & ) const;
private:
    double ct, etl;
    // cumulative time, weighted cumulative time, estimated time left, historic weight
    int n, N; // steps taken, total amount of steps
    time_t tic;
    time_t toc;
    // statics...
    static const int secperday = 86400;
    static const int secperhour = 3600;
    static const int secperminute = 60;
};
 
std::ostream & operator<<(std::ostream & , const EtaEstimator & );

#endif
