#ifndef DYN_MOD_H
#define DYN_MOD_H

#include "coordTrans.hpp"
#include <Eigen/Dense>
#include <functional>

using namespace std;
using namespace Eigen;

// State function type:  dx/dt = f(t,x,w)
typedef function<VectorXd(double, const VectorXd &, const VectorXd &)> StateFun;

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

// RKF78 integrator
VectorXd RKF78(const StateFun &dynODE,
               const VectorXd &initialValue,
               double from,
               double to,
               double dh);

class DynamicModel
{

public:
    // State function: dx/dt = f(t,x,w)
    const StateFun dynODE;

    // State dimension
    const int n;

    // Absolute & relative tolerance
    const double absTol, relTol;

    // Type alias for StateFun
    using stf = StateFun;

    // Constructor
    DynamicModel(const StateFun &dynODE_, int n_,
                 double absTol_, double relTol_);

    // Propagate state from ti to tf with noise w
    VectorXd operator()(double ti, double tf,
                        const VectorXd &xi, const VectorXd &w);
};

#endif
