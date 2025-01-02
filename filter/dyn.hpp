#ifndef DYN_MOD_H
#define DYN_MOD_H

#include "coordTrans.hpp"
#include <Eigen/Dense>
#include <functional>

using namespace std;
using namespace Eigen;

class DynamicModel
{

public:
    // State function type
    typedef function<VectorXd(double, const VectorXd &, const VectorXd &)> stf;

    // typedef function<VectorXd(const VectorXd &, double)> stg;

    // State function: dx/dt = f(t,x,w)
    const stf f;

    // // Transformation
    // const stg g;

    // State dimension
    const int n;

    // Absolute & relative tolerance
    const double abstol, reltol;

    // Constructor
    DynamicModel(const stf &f_, int n_,
                 double abstol_, double reltol_);

    // // Constructor
    // DynamicModel(const stg &g_);

    // // define my ODE system
    // void my_system(const state_type &x, state_type &dxdt, const double t);

    // Propagate state from ti to tf with noise w
    VectorXd operator()(double ti, double tf,
                        const VectorXd &xi, const VectorXd &w);

    // // Transfer the coordinate from ECI to MEE
    // VectorXd operator()(VectorXd &);

private:
    // Workspace
    VectorXd work;
    int iwork[5];
};

#endif
