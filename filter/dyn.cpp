#include "dyn.hpp"
#include <iostream>
#include <cmath>
#include <functional>

using namespace Eigen;
using namespace std;

VectorXd RKF78(
    const StateFun &dynODE,
    const VectorXd &initialValue,
    double from,
    double to,
    double dh)
{

    // Compute the number of steps
    size_t N = fabs(to - from) / dh;
    N += (fabs(to - from) / dh - N) > 1e-5 ? 1 : 0;
    char direction = ((to - from > 0) ? 1 : -1);

    VectorXd y = initialValue;
    MatrixXd k(y.size(), RKF78_A_TABLE.rows());

    for (size_t i = 0; i < N; ++i)
    {
        double x = direction * i * dh + from;

        // Compute k values
        for (size_t j = 0; j < RKF78_A_TABLE.rows(); ++j)
        {
            VectorXd summaryAK = VectorXd::Zero(y.size());

            for (size_t l = 0; l < j; ++l)
            {
                summaryAK += RKF78_A_TABLE(j, l) * k.col(l);
            }

            k.col(j) = dh * dynODE(x + RKF78_C_TABLE(j) * dh, y + summaryAK, VectorXd::Zero(y.size()));
        }

        // Compute the next step for y
        for (size_t j = 0; j < RKF78_B_TABLE.size(); ++j)
        {
            y += k.col(j) * RKF78_B_TABLE(j);
        }
    }
    return y; // Return the final value of y
}

DynamicModel::DynamicModel(const StateFun &dynODE_, int n_, double absTol_, double relTol_)
    : dynODE(dynODE_), n(n_), absTol(absTol_), relTol(relTol_) {}

VectorXd DynamicModel::operator()(double ti, double tf, const VectorXd &xi, const VectorXd &w)
{
    double dt = tf - ti;
    VectorXd xd = xi;

    // Use the custom RKF78 integrator
    xd = RKF78(dynODE, xi, ti, tf, dt);

    return xd + w;
}