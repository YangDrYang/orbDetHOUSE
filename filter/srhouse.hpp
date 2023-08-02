#ifndef SRHOUSE_H
#define SRHOUSE_H

#include "coordTrans.hpp"
#include "constants.hpp"
#include "house.hpp"
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

#define NO_MEASUREMENT 99999999

class SRHOUSE : public HOUSE
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

    // Constructor
    SRHOUSE(
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
    vector<Dist> distx;

    // Times
    vector<double> t;

    // Maximum time step for dynamical model
    const double dtMax;

    // Reset filter
    void reset(double t0, const Dist &distx0);

    // Save results
    void save(const string &filename, string stateType);
};

#endif
