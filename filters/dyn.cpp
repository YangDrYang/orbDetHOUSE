#include "dyn.hpp"
#include "ode.hpp"

#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <vector>

using namespace Eigen;

// Constructor
DynamicModel::DynamicModel(const stf& f_, int n_, double abstol_,
        double reltol_) : f(f_), n(n_), abstol(abstol_), reltol(reltol_),
        work(100+21*n_) {}

// Propagate state from ti to tf with noise w
Eigen::VectorXd DynamicModel::operator() (double ti, double tf,
        const Eigen::VectorXd& xi, const Eigen::VectorXd& w) {

    VectorXd x = xi;
    // typedef std::vector< double > state_type;
    if (tf > ti) {
        // auto state = [this](state_type &x, state_type &dxdt, const double t) -> void {
        //     f(t, )
        // };

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

