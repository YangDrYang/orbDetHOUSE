#ifndef MULTIPEARSON_HPP
#define MULTIPEARSON_HPP

#include <Eigen/Dense>
#include <iostream>
#include <random>
#include <cmath>
#include <vector>

using namespace std;
using namespace Eigen;

class typeIV
{
public:
    double a, r, gamma, lambda;

    uniform_real_distribution<double> unid;
    exponential_distribution<double> expd;

    typeIV(double mean, double std, double skew, double kurt);

    template <class URNG>
    double operator()(URNG &urng)
    {
        double y, R, phi;

        do
        {
            R = expd(urng);
            phi = fmod(R / gamma, M_PI);
            y = unid(urng);
        } while (y >= pow(sin(phi), r));

        return a * tan(phi - M_PI / 2) + lambda;
    }
};

class multiStatePearsonator
{
public:
    vector<typeIV> states;

    multiStatePearsonator(const VectorXd &means,
                          const VectorXd &stds,
                          const VectorXd &skews,
                          const VectorXd &kurts);

    template <class URNG>
    VectorXd operator()(URNG &urng, size_t stateIdx)
    {
        VectorXd noise = VectorXd::Zero(stateIdx);
        for (size_t i = 0; i < stateIdx; ++i)
        {
            // cout << "states[i]:\t" << states[i].a << " " << states[i].r << " " << states[i].lambda << " " << states[i].gamma << endl;
            // cout << "states[i](urng):\t" << states[i](urng) << endl;
            noise(i) = states[i](urng);
        }

        return noise;
    }
};

#endif // MULTIPEARSON_HPP
