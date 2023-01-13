#ifndef PEARSONATOR_H
#define PEARSONATOR_H

#include <random>
#include <cmath>

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

namespace Pearsonator {

    class TypeIV {

        public:

            double a, r, gamma, lambda;

            std::uniform_real_distribution<double> unid;
            std::exponential_distribution<double> expd;

            TypeIV(double mean, double std, double skew, double kurt);

            template <class URNG> double operator() (URNG& urng) {

                double y, R, phi;

                do {
                    R = expd(urng);
                    phi = fmod(R/gamma, M_PI);
                    y = unid(urng);
                } while (y >= pow(sin(phi), r));

                return a * tan(phi - M_PI/2) + lambda;

            }

    };

}

#endif
