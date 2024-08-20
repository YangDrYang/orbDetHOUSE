#include "dyn.hpp"
#include <iostream>
#include <cmath>
#include <functional>

using namespace Eigen;
using namespace std;

DynamicModel::DynamicModel(const stf &f_, int n_, double abstol_, double reltol_)
    : f(f_), n(n_), abstol(abstol_), reltol(reltol_), work(VectorXd::Zero(n_)) {}

// Reference: https://github.com/kofes/satellite-model/tree/master/src/helpers/integrals
//  Initialize RKF78_A_TABLE as an Eigen::MatrixXd
const MatrixXd RKF78_A_TABLE = (MatrixXd(13, 12) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                2.0 / 27.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                1.0 / 36.0, 1.0 / 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                1.0 / 24.0, 0.0, 1.0 / 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                5.0 / 12.0, 0.0, -25.0 / 16.0, 25.0 / 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                1.0 / 20.0, 0.0, 0.0, 1.0 / 4.0, 1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                -25.0 / 108.0, 0.0, 0.0, 125.0 / 108.0, -65.0 / 27.0, 125.0 / 54.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                31.0 / 300, 0.0, 0.0, 0.0, 61.0 / 225.0, -2.0 / 9.0, 13.0 / 900.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                2.0, 0.0, 0.0, -53.0 / 6.0, 704.0 / 45.0, -107.0 / 9.0, 67.0 / 90.0, 3.0, 0.0, 0.0, 0.0, 0.0,
                                -91.0 / 108.0, 0.0, 0.0, 23.0 / 108.0, -976.0 / 135.0, 311.0 / 54.0, -19.0 / 60.0, 17.0 / 6.0, -1.0 / 12.0, 0.0, 0.0, 0.0,
                                2383.0 / 4100, 0.0, 0.0, -341.0 / 164.0, 4496.0 / 1025.0, -301.0 / 82.0, 2133.0 / 4100.0, 45.0 / 82.0, 45.0 / 164.0, 18.0 / 41.0, 0.0, 0.0,
                                3.0 / 205.0, 0.0, 0.0, 0.0, 0.0, -6.0 / 41.0, -3.0 / 205.0, -3.0 / 41.0, 3.0 / 41.0, 6.0 / 41.0, 0.0, 0.0,
                                -1777.0 / 4100.0, 0.0, 0.0, -341.0 / 164.0, 4496.0 / 1025.0, -289.0 / 82.0, 2193.0 / 4100.0, 51.0 / 82.0, 33.0 / 164.0, 12.0 / 41.0, 0.0, 1.0)
                                   .finished();
const VectorXd RKF78_B_TABLE = (VectorXd(13) << 0.0, 0.0, 0.0, 0.0, 0.0, 34.0 / 105.0, 9.0 / 35.0, 9.0 / 35.0, 9.0 / 280.0, 9.0 / 280.0, 0, 41.0 / 840.0, 41.0 / 840.0).finished();
const VectorXd RKF78_C_TABLE = (VectorXd(13) << 0.0, 2.0 / 27.0, 1.0 / 9.0, 1.0 / 6.0, 5.0 / 12.0, 0.5, 5.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 3.0, 1.0, 0.0, 1.0).finished();

VectorXd rkf78(
    function<VectorXd(const VectorXd &, double)> f,
    const VectorXd &initialValue,
    double from,
    double to,
    double dh)
{
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

            k.col(j) = dh * f(y + summaryAK, x + RKF78_C_TABLE(j) * dh);
        }

        // Compute the next step for y
        for (size_t j = 0; j < RKF78_B_TABLE.size(); ++j)
        {
            y += k.col(j) * RKF78_B_TABLE(j);
        }
    }

    return y; // Return the final value of y
}

VectorXd DynamicModel::operator()(double ti, double tf, const VectorXd &xi, const VectorXd &w)
{
    double dt = tf - ti;
    VectorXd xd = xi;

    // Define the system function
    auto orb_dyn = [&](const VectorXd &x, double t)
    {
        VectorXd ww = VectorXd::Zero(x.size());
        return f(t, x, ww);
    };

    // Use the custom RKF78 integrator
    xd = rkf78(orb_dyn, xi, ti, tf, dt);

    return xd + w;
}