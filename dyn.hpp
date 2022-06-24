#ifndef DYN_MOD_H
#define DYN_MOD_H

#include <Eigen/Dense>
#include <functional>

class DynamicModel {

    public:

        // State function type
        typedef std::function<Eigen::VectorXd
            (double, const Eigen::VectorXd&, const Eigen::VectorXd&)> stf;

        // State function: dx/dt = f(t,x,w)
        const stf f;

        // State dimension
        const int n;

        // Absolute & relative tolerance
        const double abstol, reltol;

        // Constructor
        DynamicModel(const stf& f_, int n_,
            double abstol_, double reltol_);

        // Propagate state from ti to tf with noise w
        Eigen::VectorXd operator() (double ti, double tf,
            const Eigen::VectorXd& xi, const Eigen::VectorXd& w);

    private:

        // Workspace
        Eigen::VectorXd work;
        int iwork[5];

};

#endif
