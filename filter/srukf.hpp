#ifndef SRUKF_H
#define SRUKF_H

#include "coordTrans.hpp"
#include "constants.hpp"
#include <Eigen/Dense>
#include <functional>
#include "ukf.hpp"
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

class SRUKF : public UKF
{

public:
    // Dynamical model type
    typedef function<VectorXd(double, double, const VectorXd &, const VectorXd &)> dyn_model;

    // Measurement model type
    typedef function<VectorXd(double, const VectorXd &)> meas_model;

    // Dimensions of state, measurement, & process noise
    const int nx, nz, nw;

    // Indicates whether process noise is additive
    const bool addw;

    // Dynamical model: x(tf) = f(ti, tf, x(ti), w)
    const dyn_model f;

    // Maximum time step for dynamical model
    const double dtMax;

    // Measurement model: z = h(t, x) + v
    const meas_model h;

    // Square root of process & measurement noise covariance
    const MatrixXd Sww, Snn;

    // Standardized sigma points & weights
    const MatrixXd Sp, Su;
    const VectorXd wp, wu;

    // Number of sigma points in prediction & update step
    const int nsp, nsu;

    // Constructor
    SRUKF(
        const dyn_model &f_,
        const meas_model &h_,
        bool addw_,
        double t0,
        double dtMax_,
        const VectorXd &xm0,
        const MatrixXd &Pxx0,
        const MatrixXd &Pww_,
        const MatrixXd &Pnn_,
        sig_type stype,
        double k);

    // Prediction step
    void predict(double tp);

    // Update step with one measurement
    void update(const VectorXd &z);

    // Run filter for sequence of measurements
    void run(const VectorXd &tz, const MatrixXd &Z);

    // Generate standardized sigma points
    static MatrixXd sigmaSt(sig_type stype, int n, double k);

    // Generate sigma point weights
    static VectorXd sigmaWt(sig_type stype, int n, double k);

    // Times
    vector<double> t;

    // State estimates
    vector<VectorXd> xest;

    // State estimate covariances
    vector<MatrixXd> Pxx;

    // Square root of state estimate covariances
    vector<MatrixXd> Sxx;

    // Reset filter
    void reset(
        double t0,
        const VectorXd &xm0,
        const MatrixXd &Sxx0);

    void save(const string &filename, string stateType);

    // Directory for CUT files
    static string cut_dir;
};

MatrixXd cholupdate(MatrixXd matL, MatrixXd matW, double alpha);
#endif
