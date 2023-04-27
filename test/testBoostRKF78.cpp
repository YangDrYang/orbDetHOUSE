#include <iostream>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>

// typedef std::vector<double> state_type;

// void lorenz(const state_type &x, state_type &dxdt, double t)
// {
//     const double sigma = 10.0;
//     const double R = 28.0;
//     const double b = 8.0 / 3.0;

//     dxdt[0] = sigma * (x[1] - x[0]);
//     dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
//     dxdt[2] = -b * x[2] + x[0] * x[1];
// }

// void write_lorenz(const state_type &x, const double t)
// {
//     cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << endl;
// }

// int main()
// {
//     state_type x(3);
//     x[0] = 1.0;
//     x[1] = 0.0;
//     x[2] = 0.0;

//     double t = 0.0;
//     double dt = 0.1;

//     auto stepper = make_controlled<runge_kutta_fehlberg78<state_type>>(1E-10, 1E-10);

//     while (t < 100.0)
//     {
//         stepper.try_step(lorenz, x, t, dt);
//         write_lorenz(x, t);
//         t += dt;
//     }

//     return 0;
// }

#include <fstream>
#include <vector>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;
using namespace boost::numeric::odeint;

int main(int argc, char *argv[])
{
    Eigen::VectorXd v;

    v.resize(2);

    typedef VectorXd state_type;

    const double gam = 0.15;

    v(0) = 1.0;
    v(1) = 1.1;

    auto harmonic_oscillator = [&](const state_type &x, state_type &dxdt, const double t)
    {
        dxdt[0] = x[1];
        dxdt[1] = -x[0] - gam * x[1];
    };
    auto printer = [&](const state_type &x, const double t)
    {
        cout << "time: " << t << " state: " << x << endl;
    };
    integrate(harmonic_oscillator, v, 0.0, 10.0, 0.01, printer);

    return 0;
}