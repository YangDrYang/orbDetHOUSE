#include "dyn.hpp"
#include "ode.hpp"

#include <iostream>

using namespace Eigen;

// Constructor
DynamicModel::DynamicModel(const stf& f_, int n_, double abstol_,
        double reltol_) : f(f_), n(n_), abstol(abstol_), reltol(reltol_),
        work(100+21*n_) {}

// Propagate state from ti to tf with noise w
Eigen::VectorXd DynamicModel::operator() (double ti, double tf,
        const Eigen::VectorXd& xi, const Eigen::VectorXd& w) {

    VectorXd x = xi;

    if (tf > ti) {

        odefun fode = [this, &w] (double t, double* y, double* yd) -> void {
            Map<VectorXd> x(y,n), xd(yd,n);
            xd = f(t, x, w);
        };

        double t = ti;
        int flag = 1;

        ode(fode, n, x.data(), t, tf, reltol, abstol, flag, work.data(), iwork);

//        if (flag != 2)
//            std::cout << "ODE Error: " << flag << std::endl;

    }

    return x;

}

