#ifndef HOUSE_H
#define HOUSE_H

#include <Eigen/Dense>

#include <vector>

using namespace Eigen;

#define NO_MEASUREMENT 99999999

class HOUSE
{

public:
    // Dynamical model type
    typedef std::function<VectorXd(double, double,
                                   const VectorXd &, const VectorXd &)>
        dyn_model;

    // Measurement model type
    typedef std::function<VectorXd(double,
                                   const VectorXd &, const VectorXd &)>
        meas_model;

    // Distribution type
    class Dist
    {

    public:
        int n;

        VectorXd mean, skew, kurt;
        MatrixXd cov, covL;

        // Generate distribution from sigma points & weights
        Dist(const MatrixXd &X, const VectorXd &w);

        // Generate zero-mean Gaussian distribution
        Dist(const MatrixXd &S);
    };

    // Augmented state sigma point type
    class Sigma
    {

    public:
        int n_state, n_noise, n_pts;

        MatrixXd state, noise;
        VectorXd wgt;

        Sigma(const Dist &distX, const Dist &distW, double delta);
    };

    // Constructor
    HOUSE(
        const dyn_model &f_,
        const meas_model &h_,
        int nz_,
        double t0,
        double dtMax_,
        const Dist &distx0,
        const Dist &distw_,
        const Dist &distv_,
        double delta_);

    // Prediction step
    void predict(double tp);

    // Update step with one measurement
    void update(const VectorXd &z);

    // Run filter for sequence of measurements
    void run(const VectorXd &tz, const MatrixXd &Z);

    // Dynamical model
    const dyn_model f;

    // Measurement model type
    const meas_model h;

    // Dimensions of state & measurement
    const int nx, nz;

    // Dimensions of process & measurement noise
    const int nw, nv;

    // Distribution of process noise
    const Dist distw;

    // Distribution of measurement noise
    const Dist distv;

    // Minimal weight at mean
    const double delta;

    // History of state estimate distributions
    std::vector<Dist> distx;

    // Times
    std::vector<double> t;

    // Maximum time step for dynamical model
    const double dtMax;

    // Reset filter
    void reset(double t0, const Dist &distx0);

    // Save results
    void save(const std::string &filename);
};

#endif
