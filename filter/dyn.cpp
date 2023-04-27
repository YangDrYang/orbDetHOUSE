#include "dyn.hpp"
#include "ode.hpp"

#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace boost::numeric::odeint;
// Defining a shorthand for the type of the mathematical state
typedef vector<double> state_type;
// Constructor
DynamicModel::DynamicModel(const stf &f_, int n_, double abstol_,
                           double reltol_) : f(f_), n(n_), abstol(abstol_), reltol(reltol_),
                                             work(100 + 21 * n_) {}

// // Propagate state from ti to tf with noise w
// VectorXd DynamicModel::operator()(double ti, double tf, const VectorXd &xi, const VectorXd &w)
// {

//     VectorXd x = xi;
//     // typedef std::vector< double > state_type;
//     if (tf > ti)
//     {
//         // auto state = [this](state_type &x, state_type &dxdt, const double t) -> void {
//         //     f(t, )
//         // };

//         odefun fode = [this, &w](double t, double *y, double *yd) -> void
//         {
//             Map<VectorXd> x(y, n), xd(yd, n);
//             xd = f(t, x, w);
//         };

//         double t = ti;
//         int flag = 1;

//         ode(fode, n, x.data(), t, tf, reltol, abstol, flag, work.data(), iwork);

//         //        if (flag != 2)
//         //            std::cout << "ODE Error: " << flag << std::endl;
//     }

//     return x;
// }

// Propagate state from ti to tf with noise w using Boost RKF78 integrator
VectorXd DynamicModel::operator()(double ti, double tf, const VectorXd &xi, const VectorXd &w)
{

    // Error stepper, used to create the controlled stepper
    typedef runge_kutta_fehlberg78<state_type> rkf78;
    // Controlled stepper:
    // it's built on an error stepper and allows us to have the output at each
    // internally defined (refined) timestep, via integrate_adaptive call
    typedef controlled_runge_kutta<rkf78> ctrl_rkck78;
    auto my_system = [&](const state_type &x, state_type &dxdt, const double t)
    {
        VectorXd ww = VectorXd::Zero(x.size());
        VectorXd xx = VectorXd::Map(&x[0], x.size());
        VectorXd dxdt0 = f(t, xx, ww);
        dxdt = *(new vector<double>(dxdt0.data(), dxdt0.data() + dxdt0.size()));
    };
    double dt = tf - ti;
    state_type xd(xi.data(), xi.data() + xi.size());
    integrate_adaptive(ctrl_rkck78(), my_system, xd, ti, tf, dt);
    return VectorXd::Map(&xd[0], xd.size());
    // VectorXd x = xi;
    // if (tf > ti)
    // {
    //     double dt = tf - ti;
    //     int ndim = n;
    //     const state_type x0(xi);
    //     integrate_adaptive(
    //         rkf78, [&](const state_type &y, state_type &dydt, double t)
    //         {
    //             VectorXd x = y;
    //             VectorXd xd = f(t, x, w);
    //             dydt = xd; },
    //         x0, ti, tf, dt);
    // }
}
